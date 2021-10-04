#' ANCOM-BC analysis for multiple groups
#' @description This function is a wrapper of \link[ANCOMBC]{ancombc}. When
#'   having multiple groups to compare, it allows the user to automatically loop
#'   over them and computes all possible combinations. Along with
#'   \link[ANCOMBC]{ancombc} results, it appends to each table the taxonomy
#'   classification and computes the bias-adjusted abundance (expressed as
#'   logarithm), as well as the standard deviation for each group.
#'
#'
#' @param input_object_phyloseq \link[phyloseq]{phyloseq-class}. For more
#'   details, check \link[phyloseq]{distance} funtion.
#' @param grouping Character string. This is used as metadata variable for
#'   reordering factors (which allows the function to loop over groups). It is
#'   also passed to \code{group} and \link{formula} arguments in
#'   \link[ANCOMBC]{ancombc} function.
#' @param ancom.options Further arguments to be passed to
#'   \link[ANCOMBC]{ancombc}. They should be included as a list (see
#'   \code{examples})
#' @param out.unclassified Whether to filter out unclassified taxa (taxa level
#'   indicated with argument \code{tax.level}). Default is \code{FALSE}
#' @param tax.level Character string. Indicates at which level taxa should be
#'   filtered out if are unclassified.
#'
#' @return Returns a list with as many arguments as groups are contained in
#'   metadata "grouping" variable. Each element is named after the group which
#'   is compared with all others and correspond to a data frame that inlcudes
#'   \link[ANCOMBC]{ancombc} parameters and the bias-adjusted abundance.
#'   In the process, it prints the groups that are being compared (first group against others).
#'
#'
#' @export
#'
#' @examples
#'
#'ancom_location<- ancomloop(location_phyloseq,grouping="location",
#'ancom.options=list(p_adj_method ="holm"), out.unclassified=TRUE, tax.level="Genus")

#'
#'
#'
#'
#'
ancomloop <-  function (input_object_phyloseq, grouping,
                        ancom.options, out.unclassified=FALSE, tax.level){

  ########################################
  #Make grouping variable a factor
  mt <- phyloseq::sample_data(input_object_phyloseq)
  mt[[grouping]] <- as.factor(mt[[grouping]])

  #Create empty variables
  list_ancom_results <- list()
  name_out <- NULL
  #######################################
  #Loop over grouping factor levels
  for (i in seq_along(levels(mt[[grouping]]))){

    #generate reorder levels with forcats
    reord_levels <-  forcats::fct_relevel(levels(mt[[grouping]]), levels(mt[[grouping]])[i], after=0)

    #Print levels to follow the process on command line
    print(levels(reord_levels))

    #Add new reorder factor to phyloseq object
    phyloseq::sample_data(input_object_phyloseq)[[grouping]] <- factor(phyloseq::sample_data(input_object_phyloseq)[[grouping]], levels = levels(reord_levels))
    phy_data_factor=phyloseq::phyloseq(phyloseq::otu_table(input_object_phyloseq), phyloseq::tax_table(input_object_phyloseq),phyloseq::sample_data(input_object_phyloseq)) #Create factor variable for comparison
    #Use ancombc on the new phyloseq object

    if (missing(ancom.options)){

      response_out  <-  do.call(ANCOMBC::ancombc, c(list(phyloseq = phy_data_factor,formula = grouping, group=grouping)))

    }
    else {response_out  <-  do.call(ANCOMBC::ancombc, c(list(phyloseq = phy_data_factor,formula = grouping, group=grouping), ancom.options))
    }
    #Use ancombc on the new phyloseq object


    ######################PROCESS RESULTS########################################

    #Save each element as data frame
    table_ancom <- response_out$res
    #Change colnames to reflect comparison and reorder it.

    #First, save each index (it containes as many columns as comparisons)
    beta <- table_ancom$beta
    se <- table_ancom$se
    W <-  table_ancom$W
    pval <- table_ancom$p_val
    qval <- table_ancom$q_val
    diff <- table_ancom$diff_abn


    #Create name for every comparison (for example,location1 vs location 2, location 1 vs location 3)
    comparison_name <-  NULL
    for(j in 1:(length(levels(reord_levels))-1)){
      comparison_name <-  append(comparison_name,paste0(levels(reord_levels)[1],"vs",levels(reord_levels)[j+1]))
    }


    #Now, paste these comparison name with index name
    colnames(beta) = paste0(comparison_name,"_beta")
    colnames(se) = paste0(comparison_name,"_SE")
    colnames(W) = paste0(comparison_name,"_W")
    colnames(pval) = paste0(comparison_name,"_pvalue")
    colnames(qval) = paste0(comparison_name,"_padjusted")
    colnames(diff) = paste0(comparison_name,"_diff")

    #Each index table is properly labeled. Then, we will paste the first element
    #of each index, in that way, we'll get all columns for each comparison
    #together

    tabla_ancom_sorted <- NULL
    for (l in 1:length(table_ancom$beta)){

      if(l==1){
        tabla_ancom_sorted <- cbind(beta[l], se[l], W[l], pval[l], qval[l], diff[l])

      } else{
        tabla_ancom_sorted <- cbind(tabla_ancom_sorted,beta[l], se[l], W[l], pval[l], qval[l], diff[l])

      }


    }


    # BIND TAXONOMY AND CORRECTED ABUNDANCES
    taxa <- phyloseq::tax_table(input_object_phyloseq)

    #GET CORRECTED ABUNDANCE FOR ANCOM
    samp_frac <-  response_out$samp_frac
    # Replace NA with 0
    samp_frac[is.na(samp_frac)] <-  0

    # Add pesudo-count (1) to avoid taking the log of 0
    log_obs_abn  <-  log(microbiome::abundances(input_object_phyloseq) + 1)
    #Adjust the log observed abundances
    log_obs_abn_adj  <-  t(t(log_obs_abn) - samp_frac)
    t_log_abn_adj <- t(log_obs_abn_adj)


    #Calculate mean and SD
    mean <- stats::aggregate(t_log_abn_adj, by=list(mt[[grouping]]), FUN=mean)%>% tibble::column_to_rownames("Group.1") %>% t()%>%
      as.data.frame()

    #Calculate OTU SD  based on grouping factor (e.g., PLOT) and change colnames
    SD <- stats::aggregate(t_log_abn_adj, by=list(mt[[grouping]]), FUN=stats::sd)%>% tibble::column_to_rownames("Group.1")  %>% t()  %>%
      as.data.frame() %>% dplyr::rename_with(.fn= ~paste0(colnames(mean), "SD"))

    mean <- mean %>% dplyr::rename_with(.fn= ~paste0(colnames(mean),"MeanLogAbun"))
    #Merge mean abundance, SD and taxonomy.
    merged_table <- merge(taxa, mean, by=0) %>%tibble::column_to_rownames("Row.names") %>%
      merge(SD, by=0) %>% tibble::column_to_rownames("Row.names")

    #################### Merge table with ANCOM Results #########################

    table_ancom_log <- merge(merged_table, tabla_ancom_sorted, by=0) %>% tibble::column_to_rownames("Row.names")

    ##IF out.unclassified set to TRUE, filter unclassified taxa at taxonomical level set by tax.out
    if (isTRUE(out.unclassified)){
      final_ancom_table<- table_ancom_log[which(table_ancom_log[[tax.level]]!="unclassified"),]
    }else {
      final_ancom_table <- table_ancom_log

    }
    list_ancom_results[[i]] <- final_ancom_table

    #Get name of first factor which is used by ancom for comparison process. This will give name to list items
    name_out<-rbind(name_out,levels(mt[[grouping]])[i])
  }
  names(list_ancom_results) <- name_out
  return(list_ancom_results)


}
