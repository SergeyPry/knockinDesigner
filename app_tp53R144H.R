#
# This is a Shiny web application to design point mutation or codon knock-ins into 
# protein-coding and non-coding genes

library(shiny)
library(stringr)
library(Biostrings)
library(seqinr)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")),
              
  headerPanel(h1("CRISPR Knock-in Designer")),
 
  sidebarLayout(
      sidebarPanel(
        
        
        HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">Mutation site and gene data</button><p></p>'),
        
        fluidRow(
        
          column(8, textInput("Mutation", label = "Mutation name (e.g. A123C)", value = "R144H")),
          column(4, textInput("gene", label = "Gene name", value = "tp53"))
        
        ),
        
        
        fluidRow(
        
          column(12,textAreaInput("CDS", "atggcgcaaaacgacagccaagagttcgcggagctctgggagaagaatttgataagtattcagcccccaggtggtggctcttgctgggacatcattaatgatgaggagtacttgccgggatcgtttgaccccaatttttttgaaaatgtgcttgaagaacagcctcagccatccactctcccaccaacatccactgttccggagacaagcgactatcccggcgatcatggatttaggctcaggttcccgcagtctggcacagcaaaatctgtaacttgcacttattcaccggacctgaataaactcttctgtcagctggcaaaaacttgccccgttcaaatggtggtggacgttgcccctccacagggctccgtggttcgagccactgccatctataagaagtccgagcatgtggctgaagtggtccgcagatgcccccatcatgagcgaaccccggatggagataacttggcgcctgctggtcatttgataagagtggagggcaatcagcgagcaaattacagggaagataacatcactttaaggcatagtgtttttgtcccatatgaagcaccacagcttggtgctgaatggacaactgtgctactaaactacatgtgcaatagcagctgcatgggggggatgaaccgcaggcccatcctcacaatcatcactctggagactcaggaaggtcagttgctgggccggaggtcttttgaggtgcgtgtgtgtgcatgtccaggcagagacaggaaaactgaggagagcaacttcaagaaagaccaagagaccaaaaccatggccaaaaccaccactgggaccaaacgtagtttggtgaaagaatcttcttcagctacattacgacctgaggggagcaaaaaggccaagggctccagcagcgatgaggagatctttaccctgcaggtgaggggcagggagcgttatgaaattttaaagaaattgaacgacagtctggagttaagtgatgtggtgcctgcctcagatgctgaaaagtatcgtcagaaattcatgacaaaaaacaaaaaagagaatcgtgaatcatctgagcccaaacagggaaagaagctgatggtgaaggacgaaggaagaagcgactctgattaa", label ="Coding DNA sequence", height = "100px")),
          column(12,textAreaInput("exon", "tattcaccggacctgaataaactcttctgtcagctggcaaaaacttgccccgttcaaatggtggtggacgttgcccctccacagggctccgtggttcgagccactgccatctataagaagtccgagcatgtggctgaagtggtccgcagatgcccccatcatgagcgaaccccggatggagata", label ="Mutation site exon sequence",  height = "50px")),
          column(12,textAreaInput("intron5", "tatctctttaaagcaccagtactaaggaataccccagtcaataatatctcatatttaatctgctttctcattaattttagtcatgattcttacattaacttgtttagtttctaccgatactaataaaaactgcttggatggattattgaacttttttttttaagtgctaaactataacaactgggtgaaacttatttttttgtaattgcag", label ="5' intron fragment (>100 bp)", height = "50px")),
          column(12,textAreaInput("intron3", "gtacagacatttttttttccatatccattcttgcatcattctaggcctgcactattaattgattttaaaccaaaatgacgatttgaaaaggtgtgttttttttttgttgttttttgccagaaactgtgattattttgtttattactatggcgagggaggcaagtgtgtgtaattaaaacgatccaactaatgttagttaaaaagc", label ="3' intron fragment (>100 bp)", height = "50px"))          
        
        ),
        
        fluidRow(
          column(6,textInput("forw_primer", "ttgcagtattcaccggacct", label ="Forward Primer")),
          column(6,textInput("rev_primer", "gcctccctcgccatagtaat", label ="Reverse Primer"))
        ),
        
                
        
        
        HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">Guide RNA parameters</button><p></p>'),
        
        fluidRow(
          column(6, textInput("sgRNA_seq", label = "sgRNA sequence", 
                              value = "atccggggttcgctcatgat")),
          column(6,
                 selectInput("oriented", label = "sgRNA orientation", 
                             choices = list("sense" = "sense", "antisense" = "anti"), selected = "anti")
        )), 
        
        fluidRow(
          column(10, selectInput("PAM", label = "Cas9 type and PAM sequence", 
                                 choices = list("Streptococcus pyogenes-NGG" ="NGG",
                                                "Streptococcus pyogenes-NRG" ="NRG",
                                                "S.pyogenes-VQR: NGA" ="NGA",
                                                "S.pyogenes-VRER: NGCG" ="NGCG",
                                                "Staphylococcus aureus: NNGRRT" = "NNGRRT"), selected = 1))),                                 
        

        
      HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">Oligo options</button><p></p>'),
        
      fluidRow(
        
      column( 6, sliderInput("leftArmLength", "The length of left arm:", min = 30, max = 100, value = 30, step = 5)),
      column( 6, sliderInput("rightArmLength", "The length of right arm:", min = 30, max = 100, value = 30, step = 5))),
      
      fluidRow(
        
        column(6,
               selectInput("orientedOligo", label = "Oligo orientation", 
                           choices = list("sense" = "sense", "antisense" = "anti"), selected = "anti"))
        
      ),
      
      actionButton("run", "Submit")
      ), # end of sidebarPanel
      
    
      mainPanel(
        
        # UI output
        uiOutput('strategy'),
        
        uiOutput('resultsHeader'),
        
        uiOutput('oligos')       
        
      )
    )
 ))


# Define server logic to design oligos for point mutation knock-in
server <- function(input, output) {
  
  # small function to return codon differences
  getCoordinatesDiffs <- function(codon1, codon2, firstPos){
    coordinates <- c()
    
    for(i in 1:3){
      if(substr(codon1, i, i) != substr(codon2, i, i)){
        coordinates <- c(coordinates, firstPos + i - 1)
      }
    }
    return(coordinates)
  }
  
  
  # define reverse genetic code
  REV_GENETIC_CODE = list()
  
  for(codon in names(GENETIC_CODE)){
    REV_GENETIC_CODE[[as.character(GENETIC_CODE[codon])]] <- c(REV_GENETIC_CODE[[as.character(GENETIC_CODE[codon])]], codon)
  }
  
  # knockinDesign function will be used to calculate all the mutations needed to be present in oligos
  
  
  # this function will calculate the positions of all relevant items in the local genomic string
  # before mutations are engineered
  # the following need to be indicated:
  
  # 0. sequence
  # 1. Codon
  # 2. sgRNA spacer
  # 3. PAM
  # 4. Primer sites
  
  strategyCoords <- reactive({
    
    # define the output structure
    coords <- list()
    
    
    # get the coordinates of the target codon in CDS
    mutString <- input$Mutation 
    codonNum = as.integer(substr(mutString, 2, nchar(mutString)-1))
    codonCDS_pos <- c(3*codonNum-2, 3*codonNum)
    
    # determine the coordinates of the target exon inside the CDS
    # collect the input sequences
    # in the future, add validation code to make sure that clean DNA sequence is provided
    CDS <- input$CDS
    exon <- input$exon
    intr5 <- input$intron5
    intr3 <- input$intron3
    
    # perform matching of the exon to the CDS and update the codon position within the exon
    # alignment is a more generic version of solving the problem that exons 
    # may have alternative nucleotides
    alignment <- matchPattern(exon, CDS)
    d1 <- start(alignment)
    
    exonCodonPos <- codonCDS_pos - d1+1
    
    ## generate the local genomic string
    
    # update the codon coordinates
    genomicCodonPos <- exonCodonPos + nchar(intr5)
    
    # make the genomic string
    genomicString <- paste(tolower(intr5), toupper(exon), tolower(intr3), sep = "")
    
    # 0. sequence
    coords[["sequence"]] = genomicString
    coords[["exon"]] = c(nchar(intr5) + 1, nchar(intr5) + nchar(exon))
    
    # 1. Codon
    coords[["codon"]] = genomicCodonPos
    
    ## locate the sgRNA spacer
    sgRNA <- input$sgRNA_seq
    
    if(input$oriented == "sense"){
      
      # align sgRNA and genomic region when they are in the same orientation
      align_sgRNA <- matchPattern(DNAString(sgRNA), genomicString)
      start_sgR  <- start(align_sgRNA)
      end_sgR  <- end(align_sgRNA)
      
      # define the PAM coordinates
      pam_coords <- c(end_sgR +1, end_sgR + 3)
      
    } else{
      sgRNA_DS <- DNAString(sgRNA)
      sgRNA_rc <- reverseComplement(sgRNA_DS)
      
      # align sgRNA and genomic region when they are in the opposite orientations
      align_sgRNA <- matchPattern(sgRNA_rc, genomicString)
      start_sgR  <- start(align_sgRNA)
      end_sgR  <- end(align_sgRNA)
      
      # define the PAM coordinates
      pam_coords <- c(start_sgR - 3, start_sgR-1 )
      
    }
    
    # 2. sgRNA spacer
    coords[["sgRNA"]] = c(start_sgR, end_sgR)
    
    # 3. PAM
    coords[["PAM"]] = pam_coords
    
    ## primers
    forw_primer = input$forw_primer
    rev_primer = input$rev_primer  
    
    # map forward primer to the genomic string
    # toupper function is used 
    align_for_primer <- matchPattern(DNAString(toupper(forw_primer)), toupper(genomicString))
    start_for  <- start(align_for_primer)
    end_for  <- end(align_for_primer)
    
    
    # map reverse primer to the genomic string
    rev_primer <- reverseComplement(DNAString(toupper(rev_primer)))
    align_rev_primer <- matchPattern(rev_primer, DNAString(toupper(genomicString)))
    start_rev  <- start(align_rev_primer)
    end_rev  <- end(align_rev_primer)
    
    
    # 4. Primer sites
    coords[["forw_primer"]] <- c(start_for, end_for)
    coords[["rev_primer"]] <- c(start_rev, end_rev)
    
    # output the result
    coords
  })  
  
  
  
  
  codonMutations <- reactive({
    
    # get the coordinates of the target codon in CDS
    mutString <- input$Mutation 
    codonNum = as.integer(substr(mutString, 2, nchar(mutString)-1))
    codonCDS_pos <- c(3*codonNum-2, 3*codonNum)
    
    # determine the coordinates of the target exon inside the CDS
    
    # collect the input sequences
    # in the future, add validation code to make sure that clean DNA sequence is provided
    CDS <- input$CDS
    exon <- input$exon
    intr5 <- input$intron5
    intr3 <- input$intron3
    
    # perform matching of the exon to the CDS and update the codon position within the exon
    # alignment is a more generic version of solving the problem that exons 
    # may have alternative nucleotides
    alignment <- matchPattern(exon, CDS)
    d1 <- start(alignment)
    
    exonCodonPos <- codonCDS_pos - d1+1
    
    ## generate the local genomic string
    # update the codon coordinates
    genomicCodonPos <- exonCodonPos + nchar(intr5)
    
    # make the genomic string
    genomicString <- paste(tolower(intr5), toupper(exon), tolower(intr3), sep = "")
    
    ## perform codon mutations
    
    # get mutant amino acid (the last character in the mutation input)
    mutAA <- substr(input$Mutation, nchar(input$Mutation), nchar(input$Mutation))
    
    # get mutant codons
    mutCodons <- REV_GENETIC_CODE[[mutAA]]
    
    
    #############################
    # perform codon substitutions
    #############################
    # in the process, we need to generate a list of the following structure:
    
    # "ID" : number => list(
    #   "site_assay" : local_genomic_string,
    #   "new_codon": 3-letter sequence,
    #   "AA_mutation": input_mutation,
    #   "codon_coords": c(a, b) # start and end of codon in the local_genomic_string,
    #   "codon_diffs_coords": numeric or vector
    
    # set up basic output
    result = list()
    
    # start a new ID variable to store the number of the element in a list
    ID <- 1
    
    for(new_codon in mutCodons){
      # replace the previous codon at its position with a new codon 
      mutSiteAssay <- paste(substr(genomicString, 1, genomicCodonPos[1]-1),  
                            toupper(new_codon), 
                            substr(genomicString, genomicCodonPos[2]+1, nchar(genomicString)),
                            sep = "")      
      # populate the list with relevant information
      result[[ID]] <- list()
      result[[ID]][["site_assay"]] <- mutSiteAssay
      result[[ID]][["new_codon"]] <- new_codon
      result[[ID]][["AA_mutation"]] <- mutAA
      
      result[[ID]][["codon_coords"]] <- genomicCodonPos
      
      # extract the wilt-type codon from the sequence and compare it to the
      # mutant codon to extract differences
      wt_codon <- substr(genomicString, genomicCodonPos[1], genomicCodonPos[2])
      
      result[[ID]][["codon_diffs_coords"]] <- getCoordinatesDiffs(wt_codon, new_codon, genomicCodonPos[1])
      
      # update the ID that serves as the first-level key for the list
      ID = ID + 1
      
    }
    
    # output the list
    result
  })
  
  
  
  ########################################
  # PAM mutations code
  ########################################
  
  PAM_mutations <- reactive({
    
    # get  the codonMutations list input
    codons_muts <- codonMutations()
    
    # get the coordinates of sgRNA and PAM which would be 
    # common to all mutations
    coords <-strategyCoords()
    
    
    # get exon and target codon coordinates
    # in order to get a complete set of codons inside the target exon
    codonPos <- coords[["codon"]]
    exonPos <- coords[["exon"]]
    exonStart <- exonPos[1]
    exonEnd <- exonPos[2]
    
    # start a list for storing coordinates of all codons in the exon 
    codon_id = 1
    all_codons <- list()
    
    # traverse exon in a forward direction 
    nextCodonStart = codonPos[2] + 1
    nextCodonEnd = codonPos[2] + 3
    
    while(nextCodonEnd <= exonEnd ){
      all_codons[[codon_id]] <- c(nextCodonStart, nextCodonEnd)
      
      # iterate all items
      nextCodonStart = nextCodonStart + 3
      nextCodonEnd = nextCodonEnd + 3
      codon_id = codon_id + 1
      
    }
    
    # traverse exon in a backward direction 
    nextCodonStart = codonPos[1] - 3
    nextCodonEnd = codonPos[1] - 1
    
    while(nextCodonStart >= exonStart ){
      all_codons[[codon_id]] <- c(nextCodonStart, nextCodonEnd)
      
      # iterate all items
      nextCodonStart = nextCodonStart - 3
      nextCodonEnd = nextCodonEnd -3
      codon_id = codon_id + 1
      
    }
    
    # get PAM coordinates and filter the codons to select those that overlap the PAM
    coordsPAM <- coords[["PAM"]]
    
    selectedCodons <-list()
    
    codonID <- 1
    
    for(i in 1:length(all_codons)){
      # filter the codons to those that overlap PAM
      
      # check other possibilities to make sure all cases are dealt with 
      
      # check for first overlapping codon  
      if( (coordsPAM[1] <= all_codons[[i]][2]) & (coordsPAM[1] >= all_codons[[i]][1]) ){
        selectedCodons[[codonID]] <- all_codons[[i]]
        codonID <- codonID + 1
      }
      
      # check for second overlapping codon  
      if((coordsPAM[2] <= all_codons[[i]][2]) & (coordsPAM[2] >= all_codons[[i]][1])){
        selectedCodons[[codonID]] <- all_codons[[i]]
        codonID <- codonID + 1
      }
      
    }
    
    ###############################################
    # codonMutations() part testing of PAM mutations
    ################################################
    # initiate an output list
    ID <- 1
    PAM_mutated = list()
    
    # generate a match pattern based on the type of Cas9 and orientation
    if(input$oriented == "sense"){
      
      # define PAM based on the user input
      
      # "NGG"
      # "NRG"
      # "NGA"
      # "NGCG"
      # "NNGRRT"
      
      # "NGG"
      if(input$PAM == "NGG"){
        PAM_pattern = "[ACTG]GG"        
      }
      
      # "NRG"
      if(input$PAM == "NRG"){
        PAM_pattern = "[ACTG][AG]G"        
      }
      
      # "NGA"
      if(input$PAM == "NGA"){
        PAM_pattern = "[ACTG]GA"        
      }
      
      # "NGCG"
      if(input$PAM == "NGCG"){
        PAM_pattern = "[ACTG]GCG"        
      }
      
      # "NNGRRT"
      if(input$PAM == "NNGRRT"){
        PAM_pattern = "[ACTG][ACTG]G[AG][AG]T"        
      }
      
      
    } else{
      
      # define PAM based on the user input
      
      # "NGG" => "CCN"
      # "NRG" => "CYN"
      # "NGA" => "TCN"
      # "NGCG" => "CGCN"
      # "NNGRRT" => "AYYCNN"
      
      # "NGG"
      if(input$PAM == "NGG"){
        PAM_pattern = "CC[ACTG]"        
      }
      
      # "NRG"
      if(input$PAM == "NRG"){
        PAM_pattern = "C[TC][ACTG]"        
      }
      
      # "NGA"
      if(input$PAM == "NGA"){
        PAM_pattern = "TC[ACTG]"        
      }
      
      # "NGCG"
      if(input$PAM == "NGCG"){
        PAM_pattern = "CGC[ACTG]"        
      }
      
      # "NNGRRT"
      if(input$PAM == "NNGRRT"){
        PAM_pattern = "A[TC][TC]C[ACTG][ACTG]"        
      }
      
      
    } # end of PAM sequence definition
    
    #################################
    # Generation of mutant versions
    
    # set up a vector for mutated site assays
    site_assays_so_far <- c()
    
    
    # iteration over all mutated genomicString versions
    for(item in codons_muts){
      
      # initialize a flag for PAM mutations
      PAM_muts_flag <- FALSE
      
      # extract the PAM string from the mutated DNA string and test it by regular expressio
      mutSiteAssay <- item[["site_assay"]]
      curPAM <- substr(mutSiteAssay, coordsPAM[1], coordsPAM[2])
      
      # perform testing whether the pattern matches the PAM string
      if(str_detect(curPAM, PAM_pattern) ){ # PAM was not affected by codon mutation
        
        # use the selectedCodons that overlap the PAM to perform replacements
        
        # selCodon contains coordinates, not strings of codon 
        for(selCodon in selectedCodons ){
          
          # make a list of codons for which we can replace the current codon
          # start a loop for all of them
          
          # get current codon as a sequence
          selCodonSeq = substr(mutSiteAssay, selCodon[1], selCodon[2])
          
          # get all possible codons for the encoded amino acid
          aa_codons = REV_GENETIC_CODE[[ GENETIC_CODE[[selCodonSeq]] ]]
          
          # get codons that are not identical to the current codon
          codons_not_same = aa_codons[aa_codons != selCodonSeq]
          
          # iterate over these non-identical codons
          for(codon_nonID in codons_not_same){
            
            # generate a new string with a replacement codon
            mutSiteAssay_new <- paste(substr(mutSiteAssay, 1, selCodon[1]- 1),codon_nonID,
                                      substr(mutSiteAssay, selCodon[2] + 1, nchar(mutSiteAssay)), 
                                      sep = "")
            
            # select the PAM region by its coordinates
            newPAM <- substr(mutSiteAssay_new, coordsPAM[1], coordsPAM[2])
            
            # test if it matches with the pattern
            
            # if yes, do nothing, just continue to the next 
            if(str_detect(newPAM, PAM_pattern) ){
              next
              
              # if No, the mutation was successful in inactivating the PAM sequence
            }else{
              
              # test if the current site mutSiteAssay_new is already in the list
              if(mutSiteAssay_new %in% site_assays_so_far){
                next
              } else{
                
                # add the site assay to the list of such items
                site_assays_so_far <- c(site_assays_so_far, mutSiteAssay_new)
                
                # generate an entry in the output data structure  
                
                # initiate the list for this ID
                PAM_mutated[[ID]] = list()
                
                # add the new sequence
                PAM_mutated[[ID]][["site_assay"]] <- mutSiteAssay_new
                
                # copy the previous contents of the list from codon mutations 
                PAM_mutated[[ID]][["new_codon"]] <- item[["new_codon"]]
                PAM_mutated[[ID]][["AA_mutation"]] <- item[["AA_mutation"]]
                PAM_mutated[[ID]][["codon_coords"]] <- item[["codon_coords"]]
                PAM_mutated[[ID]][["codon_diffs_coords"]] <- item[["codon_diffs_coords"]]
                
                # Add the information about the PAM mutation
                PAM_mutated[[ID]][["PAM_mutant_codon"]] <- codon_nonID
                PAM_mutated[[ID]][["PAM_mut_codon_coords"]] <- selCodon
                PAM_mutated[[ID]][["PAM_mut_codon_diffs"]] <- getCoordinatesDiffs(codon_nonID,selCodonSeq, selCodon[1])
                PAM_mutated[[ID]][["sgRNA_mutations"]] <- FALSE
                
                # update the ID counter
                ID <- ID + 1
                
                # convert the flag to TRUE
                PAM_muts_flag <- TRUE
                
                
              } # end of test for duplication of site assays
              
            } # end of test whether the PAM was changed by a mutation 
          } # end non-identical codons for loop
        } # end for loop over codons overlapping with PAM
        
        
        ###############################
        # sgRNA spacer mutations
        ###############################
        
        # a flag to indicate that the sgRNA spacer was mutated
        sgRNA_mutant = FALSE
        
        # initiate empty codon vectors, which will be reassigned in the following code if 
        # some conditions are met
        
        OverlapCodon1 <- c()
        OverlapCodon2 <- c()
        
        # no replacements were successful in inactivating PAM, need to make sgRNA spacer mutations
        if(!PAM_muts_flag){
          
          # perform a search for overlapping codons based on sgRNA orientation        
          if(input$oriented == "antisense"){
            
            # iterate over all_codons inside 
            for(codon in all_codons){
              
              # a simple condition to find a codon with a minimum coordinate
              if( (codon[1] >= coords$sgRNA[1]) & (codon[1]-3 < coords$sgRNA[1]) ){
                
                OverlapCodon1 <- codon
                
                # check the second codon is within the exon
                tempNextCodon <- codon + 3
                
                if( (tempNextCodon[1] >= coords$exon[1]) & (tempNextCodon[2] <= coords$exon[2]) ){
                  OverlapCodon2 <- codon + 3
                }
                
                
              }
              
            } # end of antisense spacer overlap codon for loop           
            
            
          } else{
            
            # iterate over all_codons inside 
            for(codon in all_codons){
              
              # a simple condition to find a codon with a maximum coordinate
              if( (codon[2] <= coords$sgRNA[2]) & (codon[2]+3 > coords$sgRNA[2]) ){
                
                OverlapCodon1 <- codon
                
                # check the second codon is within the exon
                tempNextCodon <- codon - 3
                
                if( (tempNextCodon[1] >= coords$exon[1]) & (tempNextCodon[2] <= coords$exon[2]) ){
                  OverlapCodon2 <- codon - 3
                }
                
              }
              
            } # end of sense spacer overlap codon for loop           
            
            
          } # end of orientation test
          
          #########################################
          # Mutating codons within an sgRNA spacer
          #########################################
          
          # both overlap codons are defined as coordinates  
          if((length(OverlapCodon1) == 2) & (length(OverlapCodon2) == 2) ){
            
            # initiate the vector for positions of differences when a codon within sgRNA spacer is mutated 
            sgRNA_mutations = c()
            
            
            # check that the first overlap codon is not the target codon
            # and mutate it, record the positions of mutations
            if(OverlapCodon1[1] != coords$codon[1]){
              
              # get current codon as a sequence
              selCodonSeq = substr(mutSiteAssay, OverlapCodon1[1], OverlapCodon1[2])
              
              # get all possible codons for the encoded amino acid
              aa_codons = REV_GENETIC_CODE[[ GENETIC_CODE[[selCodonSeq]] ]]
              
              # get codons that are not identical to the current codon
              codons_not_same = aa_codons[aa_codons != selCodonSeq]           
              
              codon_nonID <- sample(codons_not_same, 1)
              
              # generate a new string with a replacement codon
              mutSiteAssay_new <- paste(substr(mutSiteAssay, 1, OverlapCodon1[1]- 1),codon_nonID,
                                        substr(mutSiteAssay, OverlapCodon1[2] + 1, nchar(mutSiteAssay)), 
                                        sep = "")
              
              sgRNA_mutations <- getCoordinatesDiffs(codon_nonID, selCodonSeq, OverlapCodon1[1])
              
            }
            
            
            # check that the second overlap codon is not the target codon
            # and mutate it, record the positions of mutations
            if(OverlapCodon2[1] != coords$codon[1]){
              
              # get current codon as a sequence
              selCodonSeq = substr(mutSiteAssay, OverlapCodon2[1], OverlapCodon2[2])
              
              # get all possible codons for the encoded amino acid
              aa_codons = REV_GENETIC_CODE[[ GENETIC_CODE[[selCodonSeq]] ]]
              
              # get codons that are not identical to the current codon
              codons_not_same = aa_codons[aa_codons != selCodonSeq]           
              
              codon_nonID <- sample(codons_not_same, 1)
              
              # generate a new string with a replacement codon
              mutSiteAssay_new <- paste(substr(mutSiteAssay_new, 1, OverlapCodon2[1]- 1),codon_nonID,
                                        substr(mutSiteAssay_new, OverlapCodon2[2] + 1, nchar(mutSiteAssay_new)), 
                                        sep = "")
              
              sgRNA_mutations <- c(sgRNA_mutations,getCoordinatesDiffs(codon_nonID, selCodonSeq, OverlapCodon2[1]))
              
            }          
            
            # switch the flag to TRUE to prevent the code below to be executed
            sgRNA_mutant <- TRUE
            
            
            # generate a new entry in the output list
            
            # initiate the list for this ID
            PAM_mutated[[ID]] = list()
            
            # add the new sequence
            PAM_mutated[[ID]][["site_assay"]] <- mutSiteAssay_new
            
            # copy the previous contents of the list from codon mutations 
            PAM_mutated[[ID]][["new_codon"]] <- item[["new_codon"]]
            PAM_mutated[[ID]][["AA_mutation"]] <- item[["AA_mutation"]]
            PAM_mutated[[ID]][["codon_coords"]] <- item[["codon_coords"]]
            PAM_mutated[[ID]][["codon_diffs_coords"]] <- item[["codon_diffs_coords"]]
            
            # add a flag to indicate that there were no mutations that affected the PAM site
            PAM_mutated[[ID]][["PAM_mutant_codon"]] = "none"
            
            # Add the information about the sgRNA site mutations
            PAM_mutated[[ID]][["sgRNA_mutations"]] <- TRUE
            PAM_mutated[[ID]][["sgRNA_mut_codon_diffs"]] <- sgRNA_mutations
            
            # update the ID counter
            ID <- ID + 1
            
          } # end of  both overlap codon existence
          
          # if no mutations were introduced into the sgRNA, output the previous mutant oligo entry
          if(!sgRNA_mutant){
            # copy the previous entry for the mutant codon
            PAM_mutated[[ID]] = item
            
            # add a flag to indicate that there were no mutations that affected the PAM site
            PAM_mutated[[ID]][["PAM_mutant_codon"]] = "none"
            PAM_mutated[[ID]][["sgRNA_mutations"]] <- FALSE
            
            ID = ID + 1
            
          } # end of sgRNA mutation test
        } # end of PAM_muts_flag test  
        
      }else{ # PAM was mutated so the result can be output as is
        
        # copy the previous entry for the mutant codon
        PAM_mutated[[ID]] = item
        
        # add a flag to indicate that there were no mutations that affected the PAM site
        PAM_mutated[[ID]][["PAM_mutant_codon"]] = "none"
        PAM_mutated[[ID]][["sgRNA_mutations"]] <- FALSE
        
        ID = ID + 1
      }
      
    } # end of result for loop
    
    
    PAM_mutated
  })    
  
  
  # adding silent mutations to introduce restriction sites
  # gets input from PAM_mutations()
  # this function will generate the final output to visualize the oligos and knock-in design
  
  REsite_silent_mutations <- reactive({
    
  })
  
  
  
  observeEvent(input$run, {
    
    # The strategy UI function will mark the relevant parts of the sequence to provide
    # an overview of the subsequent knock-in strategy
    output$strategy <- renderUI({
      
      # obtain the data on the strategy
      # make sure the reactive code only runs when you press "Submit" button
      isolate({ coords <- strategyCoords() })
      
      # initialize the string for HTML output
      outputHTML <- '<button type="button" class="btn btn-primary" style="width: 100%; font-size: 20px">Targeting strategy outline:</button><div class="jumbotron", style="width: 100%; word-wrap:break-word; display:inline-block;">'
      
      # process the input data to add codes to the codon, sgRNA and PAM positions
      # the primer and intervening positions can be added in a regular way
      
      for_primer_pos <- coords[["forw_primer"]]
      rev_primer_pos <- coords[["rev_primer"]]
      
      sequence = coords[["sequence"]]
      
      
      # generate the vectors of all important positions
      codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
      pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
      sgRNA_pos <- coords[["sgRNA"]][1]:coords[["sgRNA"]][2]
      
      all_special_pos <- sort(unique(c(codon_pos, pam_pos, sgRNA_pos)))
      
      
      # add forward primer stuff
      outputHTML <- paste(outputHTML,
                          "<strong><font style='BACKGROUND-COLOR: #90ee90; color: black'>",
                          substr(sequence, for_primer_pos[1], for_primer_pos[2]),
                          "</font></strong>",
                          sep = "" )
      
      # add the sequence between the forward primer and the labeled positions
      outputHTML <- paste(outputHTML,substr(sequence, for_primer_pos[2] + 1, all_special_pos[1]-1), sep="")
      
      # iterate over all labeled positions and add their 
      for(i in all_special_pos[1]: all_special_pos[length(all_special_pos)]){
        
        # Positions to be labeled
        if(i %in% all_special_pos){
          
          # Codon but NOT PAM or sgRNA
          if( (i %in% codon_pos) && !(i %in% pam_pos) && !(i %in% sgRNA_pos)){
            
            # simple yellow background of letter
            outputHTML <- paste(outputHTML,"<strong><font style='BACKGROUND-COLOR: yellow'>",
                                substr(sequence, i,i), "</font></strong>", sep="")            
          }
          
          # codon and PAM
          if( (i %in% codon_pos) && (i %in% pam_pos) && !(i %in% sgRNA_pos)){
            
            # yellow background of letter + underlined text + blue text
            outputHTML <- paste(outputHTML,"<u><strong><font style='BACKGROUND-COLOR: yellow; color: #000080'>",
                                substr(sequence, i,i), "</font></strong></u>", sep="")
          }
          
          # Codon and sgRNA
          if( (i %in% codon_pos) && !(i %in% pam_pos) && (i %in% sgRNA_pos)){
            
            # simple yellow background of letter and fuchsia-colored font
            outputHTML <- paste(outputHTML,"<strong><font style='BACKGROUND-COLOR: yellow; color: fuchsia'>",
                                substr(sequence, i,i), "</font></strong>", sep="")            
          }
          
          # PAM NOT codon
          if( (i %in% pam_pos) && !(i %in% codon_pos) && !(i %in% sgRNA_pos)){
            
            # underlined + blue text
            outputHTML <- paste(outputHTML,"<u><strong><font style='color: #000080'>",
                                substr(sequence, i,i), "</font></strong></u>", sep="")
          }
          
          # sgRNA NOT codon
          if( (i %in% sgRNA_pos) && !(i %in% codon_pos) && !(i %in% pam_pos)){
            
            # fuchsia-colored font
            outputHTML <- paste(outputHTML,"<strong><font style='color: fuchsia'>",
                                substr(sequence, i,i), "</font></strong>", sep="")            
          }
          
          
          
        } else{ # position is not labeled
          outputHTML <- paste(outputHTML,substr(sequence, i,i), sep="")
        } # end of the if statements
        
      } # end of for loop to iterate over all relevant positions
      
      
      
      # add reverse primer stuff and the rest of the sequence
      outputHTML <- paste(outputHTML, 
                          substr(sequence, all_special_pos[length(all_special_pos)] + 1, rev_primer_pos[1]-1),
                          "<strong><font style='BACKGROUND-COLOR: #90ee90; color: black'>",
                          substr(sequence, rev_primer_pos[1], rev_primer_pos[2]),
                          "</font></strong>","<br/>", sep = "")
      
      # add legend and finish this container
      outputHTML <- paste(outputHTML, "<br/>", "<p style='font-size: 12pt; font-weight:bold'>Legend: ", 
                          "<strong><font style='BACKGROUND-COLOR: #90ee90; color: black'>", 
                          "primer","</font></strong>","  ", 
                          "<strong><font style='BACKGROUND-COLOR: yellow'>", "Codon",
                          "</font></strong>","  ",
                          "<u><strong><font style='color: #000080'>", "PAM sequence",
                          "</font></strong></u>","  ",
                          "<strong><font style='color: fuchsia'>", "sgRNA sequence",
                          "</font></strong>","  ",
                          "</p>",
                          "</div>", sep = "")
      
      
      
      
      
      
      HTML(outputHTML)
    }) # end of renderUI function
  }) # end of observeEvent
  
  
  # The UI function will have to run the knockinDesign function to produce the complete description for all 
  # the designs that are possible with the current input 
  # it will then take the output of the knockinDesign and present all the individual designs
  
  observeEvent(input$run, {
    
    
    output$resultsHeader <- renderUI({
      # initialize the string for HTML output
      outputHTML <- '<button type="button" class="btn btn-primary" style="width: 100%; font-size: 20px">Results of the oligo design:</button>'
      
      HTML(outputHTML)
    })  
    
    
    output$oligos <- renderUI({
      
      # obtain the output
      isolate({ outputList <- PAM_mutations() })
      
      # obtain the data on the overall strategy
      # make sure the reactive code only runs when you press "Submit" button
      isolate({ coords <- strategyCoords() })
      
      # generate the vectors of all important positions
      codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
      pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
      
      # generate coordinates for oligo start and end 
      
      
      # iterate each oligo design
      lapply(1:length(outputList), function(i) {
        
        # get sequence
        sequence <- toupper(outputList[[i]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        PAM_muts <- c()
        sgRNA_muts <- c()
        
        # codon mutations
        codon_muts <- outputList[[i]][["codon_diffs_coords"]]
        
        # get PAM mutations is they are available
        if(outputList[[i]][["PAM_mutant_codon"]] != "none"){
          
          PAM_muts <- outputList[[i]][["PAM_mut_codon_diffs"]]
          
        }
        
        # get sgRNA mutations
        if(outputList[[i]][["sgRNA_mutations"]]){
          
          sgRNA_muts <- outputList[[i]][["sgRNA_mut_codon_diffs"]] 
          
        }     
        
        # combine all mutations
        mutated <- c(codon_muts, PAM_muts, sgRNA_muts)
        
        
        
        # reverseComplement the sequence if the orientation is different
        if(input$orientedOligo == "anti"){
          
          # update the sequence
          sequence <- toString(reverseComplement(DNAString(sequence)))
          
          # update all position vectors 
          # new_pos = nchar(sequence) - pos + 1
          codon_pos <- nchar(sequence) - rev(codon_pos) + 1
          pam_pos <- nchar(sequence) - rev(pam_pos) + 1
          
          # update the mutation positions
          mutated <- nchar(sequence) - rev(mutated) + 1
          codon_muts <- nchar(sequence) - rev(codon_muts) + 1
          PAM_muts <- nchar(sequence) - rev(PAM_muts) + 1
          sgRNA_muts <- nchar(sequence) - rev(sgRNA_muts) + 1      
          
        }
        
        # organize all mutated positions into a single vector
        all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        
        # define the start and end of oligos for the purposes of correct output
        if(input$orientedOligo == "sense"){
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[1] - 4 - input$leftArmLength
          oligoEnd <- pam_pos[1] - 4 + input$rightArmLength
          
        }else{
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[2] + 3 - input$leftArmLength 
          oligoEnd <- pam_pos[2] + 3 + input$rightArmLength 
          
          
        }
        
        # consider checking whether start and end of the oligo are less and more than any labeled positions
        # in the sequence and then updating them accordingly
        # Also check that neither of the coordinates is negative or beyond the sequence length,
        # change them accordingly
        
        # TO DO: Develop a robust tabbed interface to show designs for different codons
        
        
        
        # initialize the output HTML
        outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block;'>",
                                 "<strong>Oligo design with ID #", i, "</strong>",  "<br/>", sep = ""))
        
        # add the sequence between the oligo start and just before any of the labeled positions
        outputHTML <- paste(outputHTML, substr(sequence, oligoStart, all_special_pos[1]-1), sep="")
        
        # iterate over all labeled positions and add their 
        for(i in all_special_pos[1]: all_special_pos[length(all_special_pos)]){
          
          # Positions to be labeled
          if(i %in% all_special_pos){
            ###########################################
            # Codon but NOT PAM 
            if( (i %in% codon_pos) && !(i %in% pam_pos) ){
              
              # simple yellow background of letter
              if( i %in% mutated){
                outputHTML <- paste(outputHTML,"<strong><font style='BACKGROUND-COLOR: yellow; color: red'>",
                                    substr(sequence, i,i), "</font></strong>", sep="")            
                
              }else{
                outputHTML <- paste(outputHTML,"<strong><font style='BACKGROUND-COLOR: yellow'>",
                                    substr(sequence, i,i), "</font></strong>", sep="")            
              }
              
              
            }
            
            ###########################################
            # codon and PAM
            if( (i %in% codon_pos) && (i %in% pam_pos) ){
              
              
              # simple yellow background of letter
              if( i %in% mutated){ # the position is mutated
                
                # yellow background of letter + underlined text + RED text
                outputHTML <- paste(outputHTML,"<u><strong><font style='BACKGROUND-COLOR: yellow; color: red'>",
                                    substr(sequence, i,i), "</font></strong></u>", sep="")
                
              }else{
                
                # yellow background of letter + underlined text + RED text
                outputHTML <- paste(outputHTML,"<u><strong><font style='BACKGROUND-COLOR: yellow; color: #000080'>",
                                    substr(sequence, i,i), "</font></strong></u>", sep="")
                
              }           
              
            }
            
            ###########################################
            # PAM NOT codon
            
            if( (i %in% pam_pos) && !(i %in% codon_pos) ){
              
              if( i %in% mutated){
                
                # underlined + blue text
                outputHTML <- paste(outputHTML,"<u><strong><font style='color: red'>", substr(sequence, i,i), "</font></strong></u>", sep="")
              }else{
                
                # underlined + blue text
                outputHTML <- paste(outputHTML,"<u><strong><font style='color: #000080'>", substr(sequence, i,i), "</font></strong></u>", sep="")
                
              }
              
              
            }
            
            ###########################################
            # neither Codon nor PAM position
            
            if( !(i %in% codon_pos) && !(i %in% pam_pos)){
              
              if( i %in% mutated){
                # red-colored font
                outputHTML <- paste(outputHTML,"<strong><font style='color: red'>", substr(sequence, i,i), "</font></strong>", sep="")                 
                
              }else{
                # unlabeled
                outputHTML <- paste(outputHTML, substr(sequence, i,i), sep="")   
                
              }
              
              
              
            }
            ############################################ 
            
            
          } else{ # position is not labeled
            outputHTML <- paste(outputHTML,substr(sequence, i,i), sep="")
          } # end of the if statements
          
        } # end of for loop to iterate over all relevant positions
        
        
        # add the sequence up to the end of oligo
        outputHTML <- paste(outputHTML, substr(sequence, all_special_pos[length(all_special_pos)] + 1, oligoEnd), sep = "")
        
        # add the final tags
        outputHTML <- paste(outputHTML, "<br/>", "</div>", sep = "")
        
        
        HTML(outputHTML)
        
      }) # end of lapply
      
    }) # end of renderUI
    
  }) # end of observeEvent
  
} # end of server

# Run the application 
shinyApp(ui = ui, server = server)

