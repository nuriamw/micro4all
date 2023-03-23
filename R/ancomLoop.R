#' ANCOM-BC analysis for multiple groups
#' @description This function is a wrapper of \link[ANCOMBC]{ancombc}. When
#'   having multiple groups to compare, it allows the user to automatically loop
#'   over them and computes all possible combinations. Along with
#'   \link[ANCOMBC]{ancombc} results, it appends to each table the taxonomy
#'   classification.
#' @param input_object_phyloseq \link[phyloseq]{phyloseq-class}. For more
#'   details, check \link[phyloseq]{distance} function.
#' @param grouping Character string. This is used as metadata variable for
#'   reordering factors (which allows the function to loop over groups). It is
#'   also passed to \code{group} and \link{formula} arguments in
#'   \link[ANCOMBC]{ancombc} function.
#' @param ancom.options Further arguments to be passed to
#'   \link[ANCOMBC]{ancombc}. They should be included as a list (see more information in
#'   \code{examples})
#' @param out.unclassified Whether to filter out unclassified taxa (taxa level
#'   indicated with argument \code{tax.level}). Default is \code{FALSE}
#' @param tax.level Character string. Must be a value from 'taxonomyRanks(). Indicates at which level taxa should be
#'   filtered out if are unclassified and at which taxonomical level ANCOMBC will be applied.
#'
#' @return Returns a list with as many arguments as groups are contained in
#'   metadata "grouping" variable. Each element is named after the group which
#'   is compared with all others and correspond to a data frame that includes
#'   \link[ANCOMBC]{ancombc} parameters and the bias-adjusted abundance.
#'   In the process, it prints the groups that are being compared (first group against others).
#'
#'
#' @export
#'
#' @examples
#'
#'ancom_location<- ancomloop(location_phyloseq,grouping="location",
#'ancom.options=list(p_adj_method ="holm", n_cl=4),
#'out.unclassified=TRUE, tax.level="Genus")
#'
#'
#'
#'
#'


ancomloop <-  function (input_object_phyloseq, grouping,
                               ancom.options, out.unclassified=FALSE, tax.level=NULL){

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

      response_out  <-  do.call(ANCOMBC::ancombc, c(list(phyloseq = phy_data_factor,formula = grouping, group=grouping, tax_level=tax.level)))

    }
    else {response_out  <-  do.call(ANCOMBC::ancombc, c(list(phyloseq = phy_data_factor,formula = grouping, group=grouping,tax_level=tax.level), ancom.options))
    }
    #Use ancombc on the new phyloseq object


    ######################PROCESS RESULTS########################################

    #Save each element as data frame
    table_ancom <- response_out$res
    #Change colnames to reflect comparison and reorder it.

    #First, save each index (it contains as many columns as comparisons)
    lfc <- table_ancom[[1]]
    se <- table_ancom[[2]]
    W <-  table_ancom[[3]]
    pval <- table_ancom[[4]]
    qval <- table_ancom[[5]]
    diff <- table_ancom[[6]]


    #Create name for every comparison (for example,location1 vs location 2, location 1 vs location 3)
    comparison_name <-  NULL
    for(j in 1:(length(levels(reord_levels))-1)){
      comparison_name <-  append(comparison_name,paste0(levels(reord_levels)[1],"vs",levels(reord_levels)[j+1]))
    }


    #Now, paste these comparison name with index name
    colnames(lfc) = c("Taxon","Intercept_lfc",paste0(comparison_name,"_",names(table_ancom)[1]))
    colnames(se) = c("Taxon","Intercept_se",paste0(comparison_name,"_",names(table_ancom)[2]))
    colnames(W) = c("Taxon","Intercept_W",paste0(comparison_name,"_",names(table_ancom)[3]))
    colnames(pval) = c("Taxon","Intercept_pval",paste0(comparison_name,"_",names(table_ancom)[4]))
    colnames(qval) = c("Taxon","Intercept_qval",paste0(comparison_name,"_",names(table_ancom)[5]))
    colnames(diff) = c("Taxon","Intercept_diff",paste0(comparison_name,"_",names(table_ancom)[6]))

    #Each index table is properly labeled. Then, we will paste the first element
    #of each index, in that way, we'll get all columns for each comparison
    #together

    tabla_ancom_sorted <- NULL
    for (l in 1:length(table_ancom$lfc)){

      if(l==1){
        tabla_ancom_sorted <- cbind(lfc[l], se[l], W[l], pval[l], qval[l], diff[l])

      } else{
        tabla_ancom_sorted <- cbind(tabla_ancom_sorted,lfc[l], se[l], W[l], pval[l], qval[l], diff[l])

      }


    }


    # BIND TAXONOMY AND CORRECTED ABUNDANCES
    ##Remove repeated "Taxon" columns
    tabla_ancom_sorted <- tabla_ancom_sorted[,6:ncol(tabla_ancom_sorted)]
    ##Set "Taxon" colname to tax.level
    colnames(tabla_ancom_sorted)[1] <- tax.level

    glom_phy <- phyloseq::tax_glom(input_object_phyloseq,
                                   taxrank = tax.level)
    taxa <- BiocGenerics::as.data.frame(phyloseq::tax_table(glom_phy))
    #################### Merge table with ANCOM Results #########################
    table_ancom_log <- merge(taxa, tabla_ancom_sorted, by=tax.level)

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
