# knockinDesigner shiny app

library(shiny)
library(shinyFeedback)
library(stringr)
library(Biostrings)
library(seqinr)
library(biomaRt)
library(DECIPHER)
library(TmCalculator)
library(shinyjs)


# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
  
  useShinyFeedback(), # include shinyFeedback
  
  useShinyjs(),
  
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
            tags$style(HTML(".shiny-output-error-validation {color: green;font-size: 20px;}"))),
  
  
  HTML('<div class="row"><div style="background: #325d88; width: 100%; height: 100px; text-indent: 10px; line-height: 80px; font-size: 35px; text-align: center; color: white; text-transform: none; "><img src="knockin_logo2.png" style="float:left; height="100px"; width="223px">CRISPR Knock-in Designer</div></div>'),
  
  sidebarLayout(
    sidebarPanel(
      
      
      HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">Gene mutation</button><p></p>'),
      
      fluidRow(
        column(5, textInput("gene", label = "Gene name (optional)", value = "tp53")),
        column(6, textInput("Mutation", label = "Mutation (e.g. A123C)", value = "R144H"))
      ),
      
      HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">Gene sequence data</button><p></p>'),
      
      tabsetPanel(type = "tabs",
                  tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Manual data input"), 
                           tags$p(),
                           
                           fluidRow(
                             
                             column(12,textAreaInput("CDS", "atggcgcaaaacgacagccaagagttcgcggagctctgggagaagaatttgataagtattcagcccccaggtggtggctcttgctgggacatcattaatgatgaggagtacttgccgggatcgtttgaccccaatttttttgaaaatgtgcttgaagaacagcctcagccatccactctcccaccaacatccactgttccggagacaagcgactatcccggcgatcatggatttaggctcaggttcccgcagtctggcacagcaaaatctgtaacttgcacttattcaccggacctgaataaactcttctgtcagctggcaaaaacttgccccgttcaaatggtggtggacgttgcccctccacagggctccgtggttcgagccactgccatctataagaagtccgagcatgtggctgaagtggtccgcagatgcccccatcatgagcgaaccccggatggagataacttggcgcctgctggtcatttgataagagtggagggcaatcagcgagcaaattacagggaagataacatcactttaaggcatagtgtttttgtcccatatgaagcaccacagcttggtgctgaatggacaactgtgctactaaactacatgtgcaatagcagctgcatgggggggatgaaccgcaggcccatcctcacaatcatcactctggagactcaggaaggtcagttgctgggccggaggtcttttgaggtgcgtgtgtgtgcatgtccaggcagagacaggaaaactgaggagagcaacttcaagaaagaccaagagaccaaaaccatggccaaaaccaccactgggaccaaacgtagtttggtgaaagaatcttcttcagctacattacgacctgaggggagcaaaaaggccaagggctccagcagcgatgaggagatctttaccctgcaggtgaggggcagggagcgttatgaaattttaaagaaattgaacgacagtctggagttaagtgatgtggtgcctgcctcagatgctgaaaagtatcgtcagaaattcatgacaaaaaacaaaaaagagaatcgtgaatcatctgagcccaaacagggaaagaagctgatggtgaaggacgaaggaagaagcgactctgattaa", label ="Coding DNA sequence", height = "100px")),
                             column(12,textAreaInput("exon", "tattcaccggacctgaataaactcttctgtcagctggcaaaaacttgccccgttcaaatggtggtggacgttgcccctccacagggctccgtggttcgagccactgccatctataagaagtccgagcatgtggctgaagtggtccgcagatgcccccatcatgagcgaaccccggatggagata", label ="Mutation site exon sequence",  height = "50px")),
                             column(12,textAreaInput("intron5", "tatctctttaaagcaccagtactaaggaataccccagtcaataatatctcatatttaatctgctttctcattaattttagtcatgattcttacattaacttgtttagtttctaccgatactaataaaaactgcttggatggattattgaacttttttttttaagtgctaaactataacaactgggtgaaacttatttttttgtaattgcag", label ="5' flanking fragment (>100 bp if possible)", height = "50px")),
                             column(12,textAreaInput("intron3", "gtacagacatttttttttccatatccattcttgcatcattctaggcctgcactattaattgattttaaaccaaaatgacgatttgaaaaggtgtgttttttttttgttgttttttgccagaaactgtgattattttgtttattactatggcgagggaggcaagtgtgtgtaattaaaacgatccaactaatgttagttaaaaagc", label ="3' flanking fragment (>100 bp if possible)", height = "50px"))          
                           )
                  ),
                  
                  tabPanel(p(class = "panel-title",  style="width: 100%, font-size: 14px; color: blue", "ID-based input"),
                           tags$p(),
                           
                           fluidRow(column(7, selectInput("species", "Species", 
                                                          choices = list("select your species" = "",
                                                                         "Human (Homo sapiens)" = "hsapiens_gene_ensembl", "Mouse (Mus musculus)" = "mmusculus_gene_ensembl",
                                                                         "Rat (Rattus norvegicus)" = "rnorvegicus_gene_ensembl", "Xenopus tropicalis" = "xtropicalis_gene_ensembl",
                                                                         "Medaka (Oryzias latipes)" = "olatipes_gene_ensembl", "Zebrafish (Danio rerio)" = "drerio_gene_ensembl",
                                                                         "Fruitfly (Drosophila melanogaster)" = "dmelanogaster_gene_ensembl", "Caenorhabditis elegans" = "celegans_gene_ensembl",
                                                                         "Baker's yeast (Saccharomyces cerevisiae)" = "scerevisiae_gene_ensembl", "Gorilla (Gorilla gorilla)" = "ggorilla_gene_ensembl", 
                                                                         "Gibbon  (Nomascus leucogenys)" = "nleucogenys_gene_ensembl","Chimpanzee (Pan troglodytes)" = "ptroglodytes_gene_ensembl", 
                                                                         "Pig (Sus scrofa)" = "sscrofa_gene_ensembl", "Sheep (Ovis aries)" = "oaries_gene_ensembl",
                                                                         "Cow (Bos taurus)" = "btaurus_gene_ensembl", "Rabbit (Oryctolagus cuniculus)" = "ocuniculus_gene_ensembl",
                                                                         "Dog (Canis familiaris)" = "cfamiliaris_gene_ensembl", "Cat (Felis catus)" = "fcatus_gene_ensembl",
                                                                         "Ferret (Mustela putorius furo)" = "mfuro_gene_ensembl", "Turkey (Meleagris gallopavo)" = "mgallopavo_gene_ensembl", "Chicken (Gallus gallus)" = "ggallus_gene_ensembl",
                                                                         "Duck (Anas platyrhynchos)" = "aplatyrhynchos_gene_ensembl", "Spotted gar (Lepisosteus oculatus)" = "loculatus_gene_ensembl", 
                                                                         "Tilapia (Oreochromis niloticus )" = "oniloticus_gene_ensembl","Tetraodon viridis" = "tnigroviridis_gene_ensembl", 
                                                                         "Platyfish (Xiphophorus maculatus)" = "xmaculatus_gene_ensembl","Fugu (Fugu rubripes)" = "trubripes_gene_ensembl", 
                                                                         "Cod (Gadus morhua)" = "gmorhua_gene_ensembl", "Ciona savignyi" = "csavignyi_gene_ensembl", 
                                                                         "Ciona intestinalis" = "cintestinalis_gene_ensembl"))), 
                                    column(5, textInput("transcriptID", "Ensembl Transcript ID")) )
                 )
                  
      ), 
      
      HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">PCR primers</button><p></p>'),
      
      helpText("The primers should define an amplicon around the mutation site",
               "and ideally used before to verify guide RNA activity."),
      
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
                           choices = list("sense" = "sense", "antisense" = "anti"), selected = "anti"))
      ), 
      
      fluidRow(
        column(10, selectInput("PAM", label = "CRISPR type and PAM sequence", 
                               choices = list("Streptococcus pyogenes-NGG" ="NGG",
                                              "Streptococcus pyogenes-NRG" ="NRG",
                                              "S.pyogenes-VQR: NGA" ="NGA",
                                              "S.pyogenes-VRER: NGCG" ="NGCG",
                                              "Staphylococcus aureus: NNGRRT" = "NNGRRT",
                                              "S.aureus KKH Cas9: NNNRRT" = "NNNRRT",
                                              "AsCas12a / LbCas12a: TTTV" = "TTTV",
                                              "IDT AsCpf1/Cas12a: TTTN" = "TTTN",
                                              "FnCas12a: TTN" = "TTN",
                                              "FnCas12a: YTN" = "YTN",
                                              "Mb3Cas12a: NTTN" = "NTTN",
                                              "ErCas12a: YTTN" = "YTTN"), selected = 1))),                                 
      
      
      
      HTML('<button type="button" class="btn btn-primary" style="width: 100%; font-size: 14px">Oligo options</button><p></p>'),
      
      fluidRow(
        
        column( 6, sliderInput("leftArmLength", "left arm length", min = 30, max = 100, value = 30, step = 5)),
        column( 6, sliderInput("rightArmLength", "right arm length", min = 30, max = 100, value = 30, step = 5))),
      
      fluidRow(
        
        column(6,
               selectInput("orientedOligo", label = "Oligo orientation", 
                           choices = list("sense" = "sense", "antisense" = "anti"), selected = "sense"))
        
      ),
      
     
      
      fluidRow(
        
        column(6, radioButtons("mutatePAM", "Synonymous codon mutations of PAM or sgRNA spacer?",
                               c("Yes" = "yes", "No" = "no")),
               helpText("PAM/sgRNA codon mutations are only suitable for Cas9 and not for Cas12a/Cpf1 enzymes.")
               ),
        column(6, radioButtons("REsites", "Introduce restriction enzyme sites by synonymous codon mutations?",
                               c("Yes" = "yes", "No" = "no")))
        
      ),
      
       
      
      actionButton("run", "Submit")
    ), # end of sidebarPanel      
    
    mainPanel(
      # a bit of space from the top banner
      tags$p(),
      
      
      # panels for instructions 
      tabsetPanel(id="maintabset",
                  
                  tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Instructions"), value = "Instructions", includeHTML("instructions.html")),
                  
                  tabPanel(p(class = "panel-title",style="width: 100%, font-size: 14px; color: blue", "Results"), value = "Results", 
                          
                           # commented out because of no useful purpose                           
                           # conditionalPanel(
                           #   condition = "input.run == FALSE", 
                           #   tags$p("Oligo output data will be shown here once the program runs.")
                           # ),
                           
                           conditionalPanel(
                             
                             condition = "input.run == TRUE",
                             
                             # UI output
                             uiOutput('finalOligos')
                             
                           )
                           
                           
                           
                  )
                  
      )
      
      

      
    ) # end of mainPanel
    
  )
))



# Define server logic to design oligos for point mutation knock-in
server <- function(input, output, session) {
  
  # this is the code for updating the relevant tab panel
  observeEvent(input$run, {
  if( (input$run >= 1) & (input$Mutation != '') ){
    
    updateTabsetPanel(session, "maintabset", selected = "Results")    
  }

  })

  
  #################################
  # shinyFeedback section
  #################################

  # shinyFeedback code for Mutation input
  observeEvent(input$Mutation, {
    # Mutation is not empty but it does not conform to the pattern
    if( (str_trim(input$Mutation) != "") & !str_detect(str_trim(input$Mutation), "^[a-zA-Z*]\\d{1,3}[a-zA-Z*]$") ){
      
      feedbackWarning(
        inputId = "Mutation",
        condition = !str_detect(str_trim(input$Mutation), "^[a-zA-Z*]\\d{1,5}[a-zA-Z*]$"),
        text = "Mutation must have the pattern A123C or A123*. Please correct this."
      )
      
    } else{ # Mutation conforms to the pattern
      
      feedbackSuccess(
        inputId = "Mutation",
        condition = str_detect(str_trim(input$Mutation), "^[a-zA-Z*]\\d{1,5}[a-zA-Z*]$"),
        text = " "
      )
      
      
    }
    
    # the input has been submitted but the Mutation input was empty      
    if(input$run & str_trim(input$Mutation) == ""){
      
      feedbackWarning(
        inputId = "Mutation",
        condition = str_trim(input$Mutation) == "",
        text = "Mutation name cannot be blank. Please enter a mutation."
      )
      
    }
    
    
    
    
  })
  
  # shinyFeedback code for CDS input
  observeEvent(input$CDS, {
    
    # the input has been submitted but the CDS input was empty      
    if(input$run & str_trim(input$CDS) == ""){
      
      feedbackWarning(
        inputId = "CDS",
        condition = str_trim(input$CDS) == "",
        text = "CDS cannot be blank. Please enter the coding sequence."
      )
      
    }
    
    # give a Success message once the CDS satisfies the requirements
    if(str_trim(input$CDS) != ""){
      CDS_input <- toupper(str_trim(input$CDS))
      lettersCDS <- unique(strsplit(CDS_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "CDS",
        condition = (nchar(CDS_input) %% 3 == 0) & (sum(lettersCDS %in% DNA_BASES) == length(lettersCDS)),
        text = " "
      )
      
      # 
      if( (nchar(CDS_input) %% 3 != 0) | (sum(lettersCDS %in% DNA_BASES) != length(lettersCDS))  ){
        
        feedbackWarning(
          inputId = "CDS",
          condition = (nchar(CDS_input) %% 3 != 0) | (sum(lettersCDS %in% DNA_BASES) != length(lettersCDS)),
          text = "CDS is not fully correct. Please provide correct CDS."
        )
        
      }
      
      
    }
    
  })
  
  # shinyFeedback code for exon input
  observeEvent(input$exon, {
    
    # the input has been submitted but the exon input was empty      
    if(input$run & str_trim(input$exon) == ""){
      
      feedbackWarning(
        inputId = "exon",
        condition = str_trim(input$exon) == "",
        text = "Exon sequence cannot be blank. Please enter the exon sequence."
      )
      
    }
    
    # give a Success message once an exon satisfies the requirements
    if(str_trim(input$exon) != ""){
      exon_input <- toupper(str_trim(input$exon))
      lettersExon <- unique(strsplit(exon_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "exon",
        condition = sum(lettersExon %in% DNA_BASES) == length(lettersExon),
        text = " "
      )
      
      # 
      if( sum(lettersExon %in% DNA_BASES) != length(lettersExon)  ){
        
        feedbackWarning(
          inputId = "exon",
          condition = sum(lettersExon %in% DNA_BASES) != length(lettersExon),
          text = "Exon sequence contains non-DNA characters. Please provide correct exon sequence."
        )
        
      }
      
      
    }
    
  })
  
  # shinyFeedback code for 5' flanking sequence input
  observeEvent(input$intron5, {
    
    # the input has been submitted but the exon input was empty      
    if(input$run & str_trim(input$intron5) == ""){
      
      feedbackWarning(
        inputId = "intron5",
        condition = str_trim(input$intron5) == "",
        text = "5' flanking sequence cannot be blank. Please enter the 5' flanking sequence."
      )
      
    }
    
    # give a Success message once a 5' flanking sequence satisfies the requirements
    if(str_trim(input$intron5) != ""){
      intr5_input <- toupper(str_trim(input$intron5))
      lettersIntr5 <- unique(strsplit(intr5_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "intron5",
        condition = sum(lettersIntr5 %in% DNA_BASES) == length(lettersIntr5),
        text = " "
      )
      
      # 
      if( sum(lettersIntr5 %in% DNA_BASES) != length(lettersIntr5)  ){
        
        feedbackWarning(
          inputId = "intron5",
          condition = sum(lettersIntr5 %in% DNA_BASES) != length(lettersIntr5),
          text = "Your 5' flanking sequence contains non-DNA characters. Please provide a correct 5' flanking sequence."
        )
        
      }
      
      
    }
    
  })
  
  # shinyFeedback code for 3' flanking sequence input
  observeEvent(input$intron3, {
    
    # the input has been submitted but the exon input was empty      
    if(input$run & str_trim(input$intron3) == ""){
      
      feedbackWarning(
        inputId = "intron3",
        condition = str_trim(input$intron3) == "",
        text = "3' flanking sequence cannot be blank. Please enter the 3' flanking sequence."
      )
      
    }
    
    # give a Success message once a 3' flanking sequence satisfies the requirements
    if(str_trim(input$intron3) != ""){
      intr3_input <- toupper(str_trim(input$intron3))
      lettersIntr3 <- unique(strsplit(intr3_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "intron3",
        condition = sum(lettersIntr3 %in% DNA_BASES) == length(lettersIntr3),
        text = " "
      )
      
      # 
      if( sum(lettersIntr3 %in% DNA_BASES) != length(lettersIntr3)  ){
        
        feedbackWarning(
          inputId = "intron3",
          condition = sum(lettersIntr3 %in% DNA_BASES) != length(lettersIntr3),
          text = "Your 3' flanking sequence contains non-DNA characters. Please provide a correct 3' flanking sequence."
        )
        
      }
    }
  })
  
  # shinyFeedback code for the forward primer sequence
  observeEvent(input$forw_primer, {
    
    # the input has been submitted but the forward primer input was empty      
    if(input$run & str_trim(input$forw_primer) == ""){
      
      feedbackWarning(
        inputId = "forw_primer",
        condition = str_trim(input$forw_primer) == "",
        text = "Forward primer sequence cannot be blank. Please enter the forward primer sequence."
      )
      
    }
    
    # give a Success message once a forward primer sequence satisfies the requirements
    if(str_trim(input$forw_primer) != ""){
      forw_primer_input <- toupper(str_trim(input$forw_primer))
      lettersForw <- unique(strsplit(forw_primer_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "forw_primer",
        condition = (sum(lettersForw %in% DNA_BASES) == length(lettersForw)) & nchar(forw_primer_input) >= 15,
        text = " "
      )
      
      # warning in case non-DNA characters are within the forward primer sequence
      if( sum(lettersForw %in% DNA_BASES) != length(lettersForw)  ){
        
        feedbackWarning(
          inputId = "forw_primer",
          condition = sum(lettersForw %in% DNA_BASES) != length(lettersForw),
          text = "Your forward primer contains non-DNA characters. Please provide a correct forward primer sequence."
        )
        
      }
      
      # warning in case the forward primer is too short
      if( nchar(forw_primer_input) < 15 ){
        
        feedbackWarning(
          inputId = "forw_primer",
          condition = nchar(forw_primer_input) < 15,
          text = "Make sure your forward primer is at least 15 nucleotides long."
        )
        
      }
      
      
    }
  })
  
  # shinyFeedback code for the reverse primer sequence
  observeEvent(input$rev_primer, {
    
    # the input has been submitted but the reverse primer input was empty      
    if(input$run & str_trim(input$rev_primer) == ""){
      
      feedbackWarning(
        inputId = "rev_primer",
        condition = str_trim(input$rev_primer) == "",
        text = "Reverse primer sequence cannot be blank. Please enter the reverse primer sequence."
      )
      
    }
    
    # give a Success message once a reverse primer sequence satisfies the requirements
    if(str_trim(input$rev_primer) != ""){
      rev_primer_input <- toupper(str_trim(input$rev_primer))
      lettersRev <- unique(strsplit(rev_primer_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "rev_primer",
        condition = (sum(lettersRev %in% DNA_BASES) == length(lettersRev)) & nchar(rev_primer_input) >= 15,
        text = " "
      )
      
      # warning in case non-DNA characters are within the reverse primer sequence
      if( sum(lettersRev %in% DNA_BASES) != length(lettersRev)  ){
        
        feedbackWarning(
          inputId = "rev_primer",
          condition = sum(lettersRev %in% DNA_BASES) != length(lettersRev),
          text = "Your reverse primer contains non-DNA characters. Please provide a correct reverse primer sequence."
        )
        
      }
      
      # warning in case the reverse primer is too short
      if( nchar(rev_primer_input) < 15 ){
        
        feedbackWarning(
          inputId = "rev_primer",
          condition = nchar(rev_primer_input) < 15,
          text = "Make sure your reverse primer is at least 15 nucleotides long."
        )
        
      }
      
      
    }
  })
  
  # shinyFeedback code for the sgRNA sequence
  observeEvent(input$sgRNA_seq, {
    
    # the input has been submitted but the sgRNA sequence input was empty      
    if(input$run & str_trim(input$sgRNA_seq) == ""){
      
      feedbackWarning(
        inputId = "sgRNA_seq",
        condition = str_trim(input$sgRNA_seq) == "",
        text = "sgRNA sequence cannot be blank. Please enter an sgRNA sequence."
      )
    }
    
    # give a Success message once an sgRNA sequence satisfies the requirements
    if(str_trim(input$sgRNA_seq) != ""){
      sgRNA_input <- toupper(str_trim(input$sgRNA_seq))
      lettersSgRNA <- unique(strsplit(sgRNA_input, "")[[1]])
      
      feedbackSuccess(
        inputId = "sgRNA_seq",
        condition = (sum(lettersSgRNA %in% DNA_BASES) == length(lettersSgRNA)) & nchar(sgRNA_input) >= 18,
        text = " "
      )
      
      # warning in case non-DNA characters are within the sgRNA sequence
      if( sum(lettersSgRNA %in% DNA_BASES) != length(lettersSgRNA)  ){
        
        feedbackWarning(
          inputId = "sgRNA_seq",
          condition = sum(lettersSgRNA %in% DNA_BASES) != length(lettersSgRNA),
          text = "Your sgRNA sequence contains non-DNA characters. Please provide a correct sgRNA sequence."
        )
        
      }
      
      # warning in case the sgRNA sequence is too short
      if( nchar(sgRNA_input) < 18 ){
        
        feedbackWarning(
          inputId = "sgRNA_seq",
          condition = nchar(sgRNA_input) < 18,
          text = "Make sure your sgRNA sequence is at least 18 nucleotides long."
        )
        
      }
      
      
    }
  })
  
  #################################
  # END shinyFeedback section
  #################################
  
  
  #####################################################################
  # Lists for storing codons for each amino acid and Restriction sites
  #####################################################################  
  
  # define reverse genetic code
  REV_GENETIC_CODE = list()
  
  for(codon in names(GENETIC_CODE)){
    REV_GENETIC_CODE[[as.character(GENETIC_CODE[codon])]] <- c(REV_GENETIC_CODE[[as.character(GENETIC_CODE[codon])]], codon)
  }
  
  # define a list with restriction site information
  enzymes <- read.csv("NEB_enzymes2.csv", sep = "\t")
  
  RE_SITES <- list()
  
  for(i in rownames(enzymes)){
    enz = as.character(enzymes[i,]$Enzyme)
    
    RE_SITES[[enz]] = list()
    RE_SITES[[enz]][["Sequence"]] <- enzymes[i,]$Sequence
    RE_SITES[[enz]][["RE_site"]] <- enzymes[i,]$RE_site
    RE_SITES[[enz]][["Cut_site"]] <- enzymes[i,]$Cut_site
    RE_SITES[[enz]][["Definition"]] <- enzymes[i,]$Definition
    
  }  
  #####################################################################
  # END of Lists section
  #####################################################################
  
  
  #################################
  # FUNCTIONS
  #################################
  getPAM_pattern <- function(oriented, PAM){
    
    # generate a match pattern based on the type of Cas9 and orientation
    if(oriented == "sense"){
      # define PAM based on the user input
      
      # "NGG","NRG", "NGA", "NGCG", "NNGRRT", "NGG"
      # "TTTV", "TTTN", "TTN", "YTN", "NTTN", "YTTN"
      
      # NGG
      if(PAM == "NGG"){ PAM_pattern = "[ACTG]GG" }
      # NRG
      if(PAM == "NRG"){ PAM_pattern = "[ACTG][AG]G"}
      # NGA
      if(PAM == "NGA"){PAM_pattern = "[ACTG]GA"}
      # NGCG
      if(PAM == "NGCG"){PAM_pattern = "[ACTG]GCG"}
      # NNGRRT
      if(PAM == "NNGRRT"){PAM_pattern = "[ACTG][ACTG]G[AG][AG]T"}
      # NNNRRT
      if(PAM == "NNNRRT"){PAM_pattern = "[ACTG][ACTG][ACTG][AG][AG]T"}
      # "TTTV"
      if(PAM == "TTTV"){PAM_pattern = "TTT[ACG]"}
      # TTTN
      if(PAM == "TTTN"){PAM_pattern = "TTT[ACTG]"}      
      # TTN 
      if(PAM == "TTN"){PAM_pattern = "TT[ACTG]"}      
      # YTN 
      if(PAM == "YTN"){PAM_pattern = "[TC]T[ACTG]"}      
      # NTTN
      if(PAM == "NTTN"){PAM_pattern = "[ACTG]TT[ACTG]"}
      # YTTN
      if(PAM == "YTTN"){PAM_pattern = "[CT]TT[ACTG]"}      
      
    } else{
      
      # define PAM based on the user input
      
      # "NGG" => "CCN",  "NRG" => "CYN"
      # "NGA" => "TCN",  "NGCG" => "CGCN"
      # "NNGRRT" => "AYYCNN"
      # "TTTV" => "[TGC]AAA"
      # "TTTN" => "[ATGC]AAA"
      # "TTN"  =>  "[ATGC]AA"
      # "YTN"  =>  "[ATGC]A[AG]"      
      # "NTTN" =>  "[ATGC]AA[ATGC]"
      # "YTTN" =>  "[ATGC]AA[AG]"
      
      # "NGG"
      if(PAM == "NGG"){ PAM_pattern = "CC[ACTG]"}
      # "NRG"
      if(PAM == "NRG"){ PAM_pattern = "C[TC][ACTG]"}
      # "NGA"
      if(PAM == "NGA"){ PAM_pattern = "TC[ACTG]"}
      # "NGCG"
      if(PAM == "NGCG"){ PAM_pattern = "CGC[ACTG]"}
      # "NNGRRT"
      if(PAM == "NNGRRT"){ PAM_pattern = "A[TC][TC]C[ACTG][ACTG]" }
      # "NNNRRT"
      if(PAM == "NNNRRT"){ PAM_pattern = "A[TC][TC][ACTG][ACTG][ACTG]" }
      # "TTTV"
      if(PAM == "TTTV"){ PAM_pattern = "[TGC]AAA" }
      # "TTTN"
      if(PAM == "TTTN"){ PAM_pattern = "[ATGC]AAA" }      
      # "TTN"
      if(PAM == "TTN"){ PAM_pattern = "[ATGC]AA" }
      # "YTN"
      if(PAM == "YTN"){ PAM_pattern = "[ATGC]A[AG]" }
      # "NTTN"
      if(PAM == "NTTN"){ PAM_pattern = "[ATGC]AA[ATGC]" }
      # "YTTN"
      if(PAM == "YTTN"){ PAM_pattern = "[ATGC]AA[AG]" }      
            
    } # end of PAM sequence definition
    
    PAM_pattern
  }
  
  
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
  
  
  # function to calculate which enzymes are non-cutters
  getNonCutters <- function(all_enzymes, dna){
    
    # initialize the vector to return
    non_cutters <- c()
    
    # convert dna to character string with capital letters
    dna <- toupper(dna)
    
    # generate a reverse complement of the input DNA
    rc_dna <- as.character(reverseComplement(DNAString(dna)))
    
    # iterate over all enzymes and check whether it cuts the input sequence
    for(enzyme in all_enzymes){
      site <- as.character(RE_SITES[[enzyme]]$RE_site)
      
      # str_detect is FALSE when the site is not present
      if(!str_detect(dna, site) & !str_detect(rc_dna, site) ){
        non_cutters <- c(non_cutters, enzyme)
      }
      
    }
    
    non_cutters 
  }
  
  
  # function to identify how many enzymes now manage to cut if one of the codons is mutated
  getCutters <- function(non_cutter_enzymes, dna_piece){
    
    # vector to store the enzymes that manage to cut the sequence after mutation
    cutters <- c()
    
    # ensure that DNA is a string 
    dna <- toupper(dna_piece)
    
    # generate a reverse complement of the input DNA
    rc_dna <- as.character(reverseComplement(DNAString(dna)))
    
    # iterate over all non-cutter enzymes
    for(enzyme in non_cutter_enzymes){
      
      site <- as.character(RE_SITES[[enzyme]]$RE_site)
      
      # str_detect is TRUE when the site is present in either strand
      if( str_detect(dna, site) | str_detect(rc_dna, site)){
        cutters <- c(cutters, enzyme)
      }
      
    } # end of non_cutter enzyme for loop
    
    cutters
  }  
  
  
  # test_list is the list of the vectors of
  # same structure as the vector vec
  vector_in_list <- function(test_list, vec){
    # basic function to compare identity of one vector 
    # with another
    f1 <- function(x,y) all(x==y)  
    
    # first, if the list is empty, return FALSE
    if(length(test_list) == 0){
      return(FALSE)
      
    } else { # check if there is the same vector in the list as the query vector
      # logical result if one of the vectors matches
      return(any(mapply(f1, test_list, list(vec))))
      
    }
    
  }
  
  # a function to obtain PAM coordinates in a DNA sequence based on the type
  # and coordinates of an agRNA
  getPAM_coords <- function(PAM, oriented, start_sgR, end_sgR){
    
    if(oriented == "sense"){
      
      # define the PAM coordinates
      # based on the type of PAM that the user selects
      if(PAM %in% c("NGG", "NRG", "NGA") ){
        pam_coords <- c(end_sgR +1, end_sgR + 3)       
      }
      
      if(PAM == "NGCG"){
        pam_coords <- c(end_sgR +1, end_sgR + 4)       
      }
      
      # SaCas9 - WT or KKH
      if(PAM == "NNGRRT" | PAM == "NNNRRT"){
        pam_coords <- c(end_sgR +1, end_sgR + 6)       
      }
      
      ################ Cas12a PAM coordinates ################
      if(PAM %in% c("TTTV", "TTTN", "NTTN", "YTTN") ){
        pam_coords <- c(start_sgR - 4, start_sgR - 1)       
      }
      
      if(PAM %in% c("TTN", "YTN") ){
        pam_coords <- c(start_sgR - 3, start_sgR - 1)       
      }  
      
    } else{
      
      # define the PAM coordinates
      # based on the type of PAM that the user selects
      if(PAM %in% c("NGG", "NRG", "NGA") ){
        pam_coords <- c(start_sgR - 3, start_sgR-1 )       
      }
      
      if(PAM == "NGCG"){
        pam_coords <- c(start_sgR - 4, start_sgR-1 )       
      }
      
      # SaCas9 - WT or KKH
      if(PAM == "NNGRRT" | PAM == "NNNRRT"){
        pam_coords <- c(start_sgR - 6, start_sgR-1 )       
      }
      
      ################ Cas12a PAM coordinates ################
      if(PAM %in% c("TTTV", "TTTN", "NTTN", "YTTN") ){
        pam_coords <- c(end_sgR + 1, end_sgR + 4)       
      }
      
      if(PAM %in% c("TTN", "YTN")){
        pam_coords <- c(end_sgR + 1, end_sgR + 3)       
      }
      
    } # end of if-else 
    
    pam_coords
  } # end of function
  
  # oligoHeader function
  oligoHeader <- function(gene, mutation, codon, new_codon, j){
    
    paste("<font style='font-family: tahoma; font-size: 12; color: #325d88'><strong>", gene, " ", mutation, " (", codon, " => ", new_codon, ") oligo ", j,"</strong></font>", sep='')
  }
  
  # function to format an oligo by labeling and coloring relevant sequences using HTML
  formatOligo <- function(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated){
    outputHTML <- ""
    
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
    outputHTML <- paste(outputHTML, substr(sequence, all_special_pos[length(all_special_pos)] + 1, oligoEnd), "<br/>", sep = "")
    
    
    outputHTML
  }
  
  # calculate fragment sizes for the output
  calculateFragments <- function(RE_SITES, enzyme, site_coords, site_oriented, sequence ){
    
    # calculate the Cut_site distance
    d <- RE_SITES[[enzyme]]$Cut_site
    
    # calculate left fragment size
    if(site_oriented == "sense"){
      # sense strand
      leftFragmentSize <- site_coords[1] - 1 + d      
    } else{
      
      # anti-sense strand
      leftFragmentSize <- site_coords[1] - 1 + nchar(substr(sequence, site_coords[1], site_coords[2])) - d  
    }
    
    # calculate the right fragment size
    rightFragmentSize <- nchar(sequence) - leftFragmentSize
    
    # generate an output string
    paste("&#9986; <strong><font style='color: darkblue'>", round(leftFragmentSize/1000, 3), "kb + ", round(rightFragmentSize/1000, 3),"kb</font></strong> fragments", sep="")
  }
  
  
  # seqInputs function which processes the input and stores them in 
  # a list for future use in a similar fashion as the built-in input list
  seqInputs <- reactive({
    
    #########################################################
    # validation
    #########################################################
    
    # mutation validation
    mutString <- str_trim(input$Mutation)
    
    validate(need( mutString != "", "Mutation name cannot be blank. Please enter a mutation."))    
    validate(need(str_detect(mutString, "^[:alpha:]\\d{1,5}[a-zA-Z*]$"), "Mutation must have the pattern A123C or A123*, * - a stop codon. Please correct this." ))
    
    # if OK, store codon number
    codonNum <- as.integer(substr(mutString, 2, nchar(mutString)-1))
    
    # initialize the list for storage of data
    userData <- list()
    
    
    ##########################################################
    # retrieve the data from BioMart or from manual text input
    
    # automatic retrieval of data if all required inputs have been provided
    if(input$species != "" & input$gene != "" & input$transcriptID != "" ){
      
      # get the basic input information
      dataset_input <- input$species
      trID <- str_trim(input$transcriptID)
      
      # validate the Ensembl transcript ID
      validate( need(str_detect(trID, "ENS[A-Z]{3}T[0-9]{11}"), "Please enter a valid Ensembl Transcript ID") )
      
      
      # query Ensembl Biomart using biomaRt package
      ensembl <-  useMart("ensembl", dataset= dataset_input)
      
      # obtain exon sequences
      exon_seqs <- getSequence(id = trID, type="ensembl_transcript_id", seqType="gene_exon", mart=ensembl)
      
      # getting unspliced transcript
      unspl_transcript <- getSequence(id = trID, type="ensembl_transcript_id", seqType="transcript_exon_intron", mart=ensembl) 
      
      # once both the exon sequences and the whole unspliced transcript are obtained,
      # it will be necessary to 
      
      # 1. Map the mutation to the target exon.
      
      # define the coding sequence (CDS)
      CDS_seq <- getSequence(id = trID, type="ensembl_transcript_id", seqType="coding", mart=ensembl)
      CDS_input <- toupper(CDS_seq$coding)
      
      
      # also get the full cDNA sequence 
      CDNA_seq <- getSequence(id = trID, type="ensembl_transcript_id", seqType="cdna", mart=ensembl)
      
      
      # step 1 algorithm:
      # 1. match all exons to the cDNA to get their coordinates within the cDNA
      exon_cDNA_matches <- str_locate_all(CDNA_seq$cdna, exon_seqs$gene_exon) 
      length(exon_cDNA_matches)
      
      # iterate over the whole list and store starts and ends sep
      df <- data.frame(matrix(unlist(exon_cDNA_matches), nrow=length(exon_cDNA_matches), byrow=T))
      
      exon_seqs$start <- df$X1
      exon_seqs$end <- df$X2
      
      # 2. match the cDNA up to and including the mutation codon against the cDNA, derive the 
      #    mutation codon coordinate based on this match.
      CDS_cDNA_match <- str_locate(CDNA_seq$cdna, substr(CDS_seq$coding, 1, 3*codonNum))
      codon_end_coord <- as.integer(CDS_cDNA_match[1,]["end"])
      
      
      # 3. Identify the exon that contains the mutation based on the above matches.
      target_exon_df <- subset(exon_seqs, exon_seqs$start < codon_end_coord  & exon_seqs$end >= codon_end_coord)
      exon_input <- toupper(target_exon_df$gene_exon)
      
      # 4. Match the target exon to the unspliced transcript sequence and extract 5' and 3' intron sequences.
      target_unspl_match <- str_locate(unspl_transcript$transcript_exon_intron, toupper(exon))
      
      target_exon_unspl_start <- as.integer(target_unspl_match[1,1])
      target_exon_unspl_end <- as.integer(target_unspl_match[1,2])
      
      # 5. Store all the results of this algorithm in variables and the list for retrieval.
      
      # test that there is enough sequence on the 5' side
      if(target_exon_unspl_start >= 301){
        
        intron5 <- substr(unspl_transcript$transcript_exon_intron, target_exon_unspl_start - 300, target_exon_unspl_start-1)
        
      } else{
        
        # retrieve 300 bp upstream and join this sequence with the unspliced sequence
        flank5 <- getSequence(id = trID, type="ensembl_transcript_id", seqType="gene_flank", upstream = 300, mart=ensembl)
        
        # add the 5' flanking region to the 5' "intron"
        intron5 <- paste(flank5$gene_flank, substr(unspl_transcript$transcript_exon_intron, 1, target_exon_unspl_start-1), sep = "")
      }
      
      # test if there is enough sequence on the 3' side
      if(nchar(unspl_transcript$transcript_exon_intron) >= target_exon_unspl_end + 300 ){
        
        intron3 <- substr(unspl_transcript$transcript_exon_intron, target_exon_unspl_end + 1, target_exon_unspl_end + 300)
        
      }else{
        # retrieve 300 bp downstream from the transcript
        flank3 <- getSequence(id = trID, type="ensembl_transcript_id", seqType="gene_flank", downstream = 300, mart=ensembl)    
        
        unspl_transcript_flanked <- paste(unspl_transcript$transcript_exon_intron, flank3$gene_flank, sep = "")
        
        intron3 <- substr(unspl_transcript_flanked, target_exon_unspl_end + 1, target_exon_unspl_start + 300)
      }
      
      # store the results of intron/flanking sequence selection 
      intron5_input <- toupper(intron5)
      intron3_input <- toupper(intron3)
      
      # manual data input
    } else{
      
      # extract the sequences from input
      CDS_input <- toupper(input$CDS)
      exon_input <- toupper(input$exon)
      intron5_input <- toupper(input$intron5)
      intron3_input <- toupper(input$intron3)
      
      ## validate the input sequences
      
      # trimming the whitespace on both sides of the sequence
      CDS_input <- str_trim(CDS_input)
      # defining the unique characters in the CDS_input
      lettersCDS <- unique(strsplit(CDS_input, "")[[1]])
      
      # validation of the sequence
      validate(need( CDS_input != "", "CDS cannot be blank. Please enter a CDS"))
      validate(need( nchar(CDS_input) %% 3 == 0, "The length of CDS must be divisible by 3. Please enter a correct coding sequence."))
      validate(need(sum(lettersCDS %in% DNA_BASES) == length(lettersCDS), "Please provide a correct coding sequence"))
      
      # validation of the mutation string
      
      # extract the codon from the CDS and the corresponding amino acid
      codonCDS_pos <- c(3*codonNum-2, 3*codonNum)
      wt_codon <- substr(CDS_input, codonCDS_pos[1], codonCDS_pos[2] )
      actAA <- toupper(substr(mutString, 1, 1))
      # perform validation
      validate(need(GENETIC_CODE[[wt_codon]] == actAA, paste("The codon at position ", codonNum, " does not code for ", actAA, sep = "")))
      
      # trimming the whitespace on both sides of the sequence
      exon_input <- str_trim(exon_input)
      # defining the unique characters in the exon_input
      lettersExon <- unique(strsplit(exon_input, "")[[1]])
      
      # validation of the sequence
      validate(need( exon_input != "", "Exon sequence cannot be blank. Please enter the exon sequence"))
      validate(need(sum(lettersExon %in% DNA_BASES) == length(lettersExon), "Please provide a correct exon sequence")) 
      
      
      # trimming the whitespace on both sides of the sequence
      intron5_input <- str_trim(intron5_input)
      # defining the unique characters in the intron5_input
      lettersIntron5 <- unique(strsplit(intron5_input, "")[[1]])
      
      # validation of the sequence
      validate(need(intron5_input != "", "5' intron sequence cannot be blank. Please enter the 5' intron sequence"))
      validate(need(sum(lettersIntron5 %in% DNA_BASES) == length(lettersIntron5), "Please provide a correct 5' intron sequence")) 
      
      
      # trimming the whitespace on both sides of the sequence
      intron3_input <- str_trim(intron3_input)
      # defining the unique characters in the intron3_input
      lettersIntron3 <- unique(strsplit(intron3_input, "")[[1]])
      
      # validation of the sequence
      validate(need(intron3_input != "", "3' intron sequence cannot be blank. Please enter the 3' intron sequence"))
      validate(need(sum(lettersIntron3 %in% DNA_BASES) == length(lettersIntron3), "Please provide a correct 3' intron sequence")) 
      
      
    }# end of else for manual data input
    
    # validation of the exon being inside the CDS
    validate( need(str_detect(CDS_input, exon_input), "Please make sure your input exon matches the coding sequence.") )
    
    # generate a genomicString variable
    genomicString <- paste(intron5_input, exon_input, intron3_input, sep="")
    
    # retrieve forward primer and reverse primers and validate them
    forw_primer <- toupper(str_trim(input$forw_primer))
    rev_primer <- toupper(str_trim(input$rev_primer))
    
    # empty primer fields
    validate(need( forw_primer != "", "Forward primer sequence cannot be blank. Please enter a forward primer sequence"))
    validate(need( rev_primer != "", "Reverse primer sequence cannot be blank. Please enter a reverse primer sequence"))
    
    # matching test for the forward primer
    validate( need(str_detect(genomicString, forw_primer), "Please make sure your forward primer matches precisely either your flanking sequence or exon.") )
    
    # make a reverse complement of the reverse primer
    rc_rev_primer <- as.character(reverseComplement(DNAString(rev_primer)))
    
    # matching test for the reverse primer
    validate( need(str_detect(genomicString, rc_rev_primer), "Please make sure your reverse primer matches precisely either your exon or flanking sequence.") )
    
    # retrieve the input sgRNA and validate it with respect to the genomic string
    # this would also validate the orientation of the sgRNA
    
    
    sgRNA <- toupper(str_trim(input$sgRNA_seq))
    validate(need(sgRNA != "", "sgRNA sequence cannot be blank. Please enter the a sgRNA sequence"))
    
    if(input$oriented == "sense"){
      validate( need(str_detect(genomicString, sgRNA), "Please make sure your sgRNA matches your sequence precisely and the orientation is correct.") )
    } else{
      
      rc_sgRNA <- as.character(reverseComplement(DNAString(sgRNA)))
      
      validate( need(str_detect(genomicString, rc_sgRNA), "Please make sure your sgRNA matches your sequence precisely and the orientation is correct.") )
      
    }
    
    
    # validate the PAM sequence inside the genomicString
    
    # 1. get coordinates for sgRNA
    
    # code to define the PAM to be used in the website
    if(input$oriented == "sense"){
      
      # align sgRNA and genomic region when they are in the same orientation
      align_sgRNA <- matchPattern(DNAString(sgRNA), genomicString)
      start_sgR  <- start(align_sgRNA)
      end_sgR  <- end(align_sgRNA)
      
      # reverse sgRNA orientation            
    } else{
      sgRNA_DS <- DNAString(sgRNA)
      sgRNA_rc <- reverseComplement(sgRNA_DS)
      
      # align sgRNA and genomic region when they are in the opposite orientations
      align_sgRNA <- matchPattern(sgRNA_rc, genomicString)
      start_sgR  <- start(align_sgRNA)
      end_sgR  <- end(align_sgRNA)
    
    } # end of PAM coordinates definition 
    
    # 2. get PAM coordinates using a function
    pam_coords <- getPAM_coords(input$PAM, input$oriented, start_sgR, end_sgR)
    
    
    # 3. extract the PAM sequence based on the above, orientation and type of CRISPR system
    PAM_seq <- substr(genomicString, pam_coords[1], pam_coords[2])
    
    # generate a regular expression for PAM pattern  
    PAM_pattern <- getPAM_pattern(input$oriented, input$PAM)  
    
    # 4. perform validation based on regular expression of extracted PAM sequences with the corresponding PAM pattern
    validate( need(str_detect(PAM_seq, PAM_pattern), "Please make sure your CRISPR system is selected correctly so the PAM matches the sequence.") )
    
    #########################################################
    # END of validation section
    #########################################################
    
    
    #############################
    # OUTPUT
    #############################
    
    # the idea here is that the sequence inputs are processed or retrieved
    # everything is validated here since this function will be executed early on
    # the rest of the inputs will be accessed from the regular input list since 
    # repeated execution of the whole function will be problematic
    
    
    # store the sequence inputs pasted by the user or retrieved from the databases
    
    userData$CDS <- CDS_input
    userData$exon_sequence <- exon_input
    userData$intron5 <- intron5_input
    userData$intron3 <- intron3_input
    
    
    # output the list
    userData
  })
  
  
  # strategyCoords will calculate the positions of all relevant items in the local genomic string
  # before mutations are engineered. The following need to be indicated:
  
  # sequence, Codon, sgRNA spacer, PAM, Primer sites
  strategyCoords <- reactive({
    
    # define the output structure
    coords <- list()
    
    seqInputs <- seqInputs()
    
    # get the coordinates of the target codon in CDS
    mutString <- str_trim(input$Mutation)
    codonNum = as.integer(substr(mutString, 2, nchar(mutString)-1))
    codonCDS_pos <- c(3*codonNum-2, 3*codonNum)
    
    # determine the coordinates of the target exon inside the CDS
    # collect the input sequences
    # in the future, add validation code to make sure that clean DNA sequence is provided
    CDS <- seqInputs$CDS
    exon <- seqInputs$exon_sequence
    intr5 <- seqInputs$intron5
    intr3 <- seqInputs$intron3
    
    # store basic data in the list
    coords[["Mutation"]] <- mutString
    coords[["CDS"]] <- CDS
    coords[["exon_sequence"]] <- exon
    coords[["intron5"]] <- intr5
    coords[["intron3"]] <- intr3
    
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
    
    # code to define the PAM to be used in the website
    if(input$oriented == "sense"){
      
      # align sgRNA and genomic region when they are in the same orientation
      align_sgRNA <- matchPattern(DNAString(sgRNA), genomicString)
      start_sgR  <- start(align_sgRNA)
      end_sgR  <- end(align_sgRNA)
      
      # reverse sgRNA orientation            
    } else{
      sgRNA_DS <- DNAString(sgRNA)
      sgRNA_rc <- reverseComplement(sgRNA_DS)
      
      # align sgRNA and genomic region when they are in the opposite orientations
      align_sgRNA <- matchPattern(sgRNA_rc, toupper(genomicString))
      start_sgR  <- start(align_sgRNA)
      end_sgR  <- end(align_sgRNA)
      
    } # end of PAM coordinates definition 
    
    
    # get PAM coordinates using a function
    pam_coords <- getPAM_coords(input$PAM, input$oriented, start_sgR, end_sgR)
    
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
  
  # code to design forward AS-PCR primers
  design_forward_primers <- function(mut_seq, wt_seq, lastCodonDifference, rev_primer){
    
    # 1. make candidate primers 18-25 nucleotide long and ending at lastCodonDifference
    lengths <- seq(18,25, by=1)
    
    # initialize the vector of primers 
    primers <- c()
    
    
    for(i in 1:8){
      primers <- c(primers, substr(mut_seq, lastCodonDifference - lengths[i] + 1, lastCodonDifference))
    }
    
    # 2. calculate the Tm of each primer in the vector 
    
    # store parameters for Tm calculation 
    ambiguous=TRUE
    userset=NULL
    variant="Primer3Plus"
    Na=0
    K=50
    Tris=10
    Mg=1.5
    dNTPs=0.2
    mismatch=TRUE
    
    # vector for storing Tms
    Temperatures <- c()
    
    for(i in 1:8){
      Temperatures[i] <- Tm_NN(primers[i], ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4",
                               tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 25,
                               dnac2 = 25, selfcomp = FALSE, Na, K, Tris, Mg, dNTPs, saltcorr = 5)
    }
    
    rev_primer_Tm <- Tm_NN(rev_primer, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4",
                           tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 25,
                           dnac2 = 25, selfcomp = FALSE, Na, K, Tris, Mg, dNTPs, saltcorr = 5)
    
    # calculate Tm differences between all candidate primers and the reverse primer
    diffs <- abs(Temperatures - rev_primer_Tm) 
    
    # identify which particular item has the minimum difference
    k <- which.min(diffs)
    
    # selected mutant primer
    selected_mut_primer <- primers[k]
    
    # selected wild-type primer
    wt_primer <- substr(wt_seq, lastCodonDifference - lengths[k] + 1, lastCodonDifference)
    
    Tm_wt_primer = Tm_NN(wt_primer, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4",
                         tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 25,
                         dnac2 = 25, selfcomp = FALSE, Na, K, Tris, Mg, dNTPs, saltcorr = 5)
    
    # return the list with the mutant and wild-type forward primers
    list(mut_forw_primer = list(sequence = selected_mut_primer, Tm = Temperatures[k]), wt_forw_primer = list(sequence = wt_primer, Tm = Tm_wt_primer),
         common_rev_primer = list(sequence = rev_primer, Tm = rev_primer_Tm))  
  }
  
  # code to design forward AS-PCR primers
  design_reverse_primers <- function(mut_seq, wt_seq, firstCodonDifference, forw_primer){
    
    # 1. make candidate primers 18-25 nucleotide long and ending at lastCodonDifference
    lengths <- seq(18,25, by=1)
    
    # initialize the vector of primers 
    primers <- c()
    
    
    for(i in 1:8){
      forw_strand_oligo <- substr(mut_seq, firstCodonDifference, firstCodonDifference + lengths[i] - 1)
      rev_strand_oligo <- toString(reverseComplement(DNAString(forw_strand_oligo)))
      primers <- c(primers, rev_strand_oligo)
    }
    
    # 2. calculate the Tm of each primer in the vector 
    
    # store parameters for Tm calculation 
    ambiguous=TRUE
    userset=NULL
    variant="Primer3Plus"
    Na=0
    K=50
    Tris=10
    Mg=1.5
    dNTPs=0.2
    mismatch=TRUE
    
    # vector for storing Tms
    Temperatures <- c()
    
    for(i in 1:8){
      Temperatures[i] <- Tm_NN(primers[i], ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4",
                               tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 25,
                               dnac2 = 25, selfcomp = FALSE, Na, K, Tris, Mg, dNTPs, saltcorr = 5)
    }
    
    forw_primer_Tm <- Tm_NN(forw_primer, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4",
                            tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 25,
                            dnac2 = 25, selfcomp = FALSE, Na, K, Tris, Mg, dNTPs, saltcorr = 5)
    
    # calculate Tm differences between all candidate primers and the reverse primer
    diffs <- abs(Temperatures - forw_primer_Tm) 
    
    # identify which particular item has the minimum difference
    k <- which.min(diffs)
    
    # selected mutant primer
    selected_mut_primer <- primers[k]
    
    # selected wild-type primer
    wt_forw_oligo <- substr(wt_seq, firstCodonDifference, firstCodonDifference + lengths[k] - 1)
    wt_primer <-  toString(reverseComplement(DNAString(wt_forw_oligo)))
    
    # return the list with the mutant and wild-type forward primers
    list(mut_rev_primer = selected_mut_primer, wt_rev_primer = wt_primer)  
    
    
    Tm_wt_primer = Tm_NN(wt_primer, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = "DNA_NN4",
                         tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1", dnac1 = 25,
                         dnac2 = 25, selfcomp = FALSE, Na, K, Tris, Mg, dNTPs, saltcorr = 5)
    
    # return the list with the mutant and wild-type forward primers
    list( common_forw_primer = list(sequence = forw_primer, Tm = forw_primer_Tm), mut_rev_primer = list(sequence = selected_mut_primer, Tm = Temperatures[k]), 
          wt_rev_primer = list(sequence = wt_primer, Tm = Tm_wt_primer))
    
  }
  
  # compile primer tables
  primerTablesOutput <- function(forward_primers, reverse_primers){
    
    outputHTML <- ""
    
    ## 3. Proceed to the output stage
    
    # list(mut_forw_primer = list(sequence = selected_mut_primer, Tm = Temperatures[k]), wt_forw_primer = list(sequence = wt_primer, Tm = Tm_wt_primer),
    #      common_rev_primer = list(sequence = rev_primer, Tm = rev_primer_Tm))
    
    
    forw_primer_table <- "<table style='width:80%'><caption style='color:blue'>Forward AS-PCR assays</caption><tr><th>Primer</th><th>Sequence</th><th>Tm</th></tr>"
    forw_primer_table <- paste(forw_primer_table,  "<tr><td>", "KI forward" ,"</td>", 
                               "<td>", forward_primers[["mut_forw_primer"]]$sequence,"</td>", 
                               "<td>", round(forward_primers[["mut_forw_primer"]]$Tm, 2), "</td>", "</tr>",
                               "<tr><td>", "Wild-type forward" ,"</td>",
                               "<td>", forward_primers[["wt_forw_primer"]]$sequence,"</td>", 
                               "<td>", round(forward_primers[["wt_forw_primer"]]$Tm, 2), "</td>", "</tr>",
                               "<tr><td>", "common reverse" ,"</td>",
                               "<td>", forward_primers[["common_rev_primer"]]$sequence,"</td>", 
                               "<td>", round(forward_primers[["common_rev_primer"]]$Tm, 2), "</td>", "</tr>",
                               "</table>", sep = "")
    
    # list( common_forw_primer = list(sequence = forw_primer, Tm = forw_primer_Tm), mut_rev_primer = list(sequence = selected_mut_primer, Tm = Temperatures[k]), 
    #       wt_rev_primer = list(sequence = wt_primer, Tm = Tm_wt_primer))
    
    rev_primer_table <- "<table style='width:80%'><caption style='color:blue'>Reverse AS-PCR assays</caption><tr><th>Primer</th><th>Sequence</th><th>Tm</th></tr>"
    rev_primer_table <- paste(rev_primer_table, "<tr><td>", "common forward" ,"</td>",
                              "<td>", reverse_primers[["common_forw_primer"]]$sequence,"</td>", 
                              "<td>", round(reverse_primers[["common_forw_primer"]]$Tm, 2), "</td></tr>",
                              "<tr><td>", "KI reverse" ,"</td>", 
                              "<td>", reverse_primers[["mut_rev_primer"]]$sequence,"</td>", 
                              "<td>", round(reverse_primers[["mut_rev_primer"]]$Tm, 2), "</td></tr>",
                              "<tr><td>", "Wild-type reverse" ,"</td>",
                              "<td>", reverse_primers[["wt_rev_primer"]]$sequence,"</td>", 
                              "<td>", round(reverse_primers[["wt_rev_primer"]]$Tm, 2), "</td></tr>",
                              "</table>", sep = "")
    
    outputHTML <- paste(outputHTML, "<table style='width: 100%'><tr>")
    outputHTML <- paste(outputHTML, "<td>", forw_primer_table, "</td>", sep = "")
    outputHTML <- paste(outputHTML, "<td>", rev_primer_table, "</td>", sep = "")
    outputHTML <- paste(outputHTML, "</tr></table>")
    
    outputHTML
  }
  
  
  ########################################
  # Codon mutations code
  ########################################
  codonMutations <- reactive({
    
    # determine the coordinates of the target exon inside the CDS
    
    # collect the input sequences
    # in the future, add validation code to make sure that clean DNA sequence is provided
    coords <- strategyCoords()
    
    # get the coordinates of the target codon in CDS
    mutString <- coords$Mutation 
    codonNum = as.integer(substr(mutString, 2, nchar(mutString)-1))
    codonCDS_pos <- c(3*codonNum-2, 3*codonNum)
    
    CDS <- coords$CDS
    exon <- coords$exon_sequence
    intr5 <- coords$intron5
    intr3 <- coords$intron3
    
    # primers
    forw_primer_pos <- coords[["forw_primer"]]
    rev_primer_pos <- coords[["rev_primer"]] 
    
    # offset all the positions - necessary to convert pattern matching coordinates to normal coordinates
    offset <- forw_primer_pos[1] - 1
    
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
    
    ############################################################
    # test if existing  codon mutations generate restriction sites,
    # add them to the output list
    ############################################################
    
    # define the wild-type site assay for checking non-cutters
    wt_seq <- coords[["sequence"]]
    
    # subset the full sequence to that amplified by primers
    wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
    
    # load restriction enzyme data
    all_enzymes <- names(RE_SITES)
    
    # compute the non-cutter enzymes
    non_cutters <- getNonCutters(all_enzymes, wt_site_assay)    
    
    ######################################
    # iteration over the result list
    ######################################
    
    for(j in 1:length(result)){
      
      #################################################
      # extract relevant parts of the current item 
      #################################################
      
      # get full genomicString sequence and subset to the sequence spanned by primers
      mut_site_assay <- substr(result[[j]][["site_assay"]], forw_primer_pos[1], rev_primer_pos[2])
      
      # make a reverse complement of the mut_site_assay
      rc_mut_site_assay <- as.character(reverseComplement(DNAString(mut_site_assay)))
      
      # get coordinates of the codon in the current list item
      codon_coords <- result[[j]][["codon_coords"]] - offset
      codon_diffs_coords <- result[[j]][["codon_diffs_coords"]] - offset
      
      # make a flag variable to keep track of finding sites that 
      # have been introduced
      SITES_EXIST_FLAG <- FALSE
      
      # subset the site assay to a small 30-nt string for testing
      cur_test_string <- substr(mut_site_assay, codon_coords[1] - 15, codon_coords[1] + 18)
      
      # run a function to get all enzymes that cut 
      cutters <- getCutters(non_cutters, cur_test_string)
      
      # iterate over enzymes that cut the test string around the codon
      if(length(cutters) > 0){
        
        # initialize the list for enzyme sites
        result[[j]][["enzymes"]] <- list()
        enzID <- 1
        
        # add all of the cutter enzymes to a list inside an existing oligo entry
        for(cutterEnzyme in cutters){
          
          # get cut site
          cut_site <- as.character(RE_SITES[[cutterEnzyme]]$RE_site)
          
          # test if the cut site is in the forward strand
          if(str_detect(toupper(mut_site_assay), cut_site)){
            
            # get coordinates
            site_coords <- c(str_locate(toupper(mut_site_assay), cut_site))    
            
            # store site orientation 
            site_oriented <- "sense"
            
          }else{
            # get coordinates
            site_coords <- c(rev(nchar(rc_mut_site_assay) - str_locate(toupper(rc_mut_site_assay), cut_site) + 1))
            
            # store site orientation 
            site_oriented <- "anti"
          }
          
          # store the restriction site information
          result[[j]][["enzymes"]][[enzID]] <- list()
          
          result[[j]][["enzymes"]][[enzID]][["RE_enzyme"]] <- cutterEnzyme
          result[[j]][["enzymes"]][[enzID]][["RE_site"]] <- RE_SITES[[cutterEnzyme]]$Sequence
          
          # storing the restriction site coordinates in the original coordinates
          result[[j]][["enzymes"]][[enzID]][["RE_site_coords"]] <- site_coords + offset
          
          # store site orientation 
          result[[j]][["enzymes"]][[enzID]][["RE_site_oriented"]] <- site_oriented
          
          # update the enzID counter
          enzID <- enzID + 1
          
        } # end of cutters for loop
        
      } # end of if statement  for testing how many enzymes cut the sequence
    
    } 
    ##################################
    # end of result iteration
    ##################################
    
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
    
    # primers
    forw_primer_pos <- coords[["forw_primer"]]
    rev_primer_pos <- coords[["rev_primer"]] 
    
    # offset all the positions - necessary to convert pattern matching coordinates to normal coordinates
    offset <- forw_primer_pos[1] - 1
    
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
    
    # generate a regular expression for PAM pattern  
    PAM_pattern <- getPAM_pattern(input$oriented, input$PAM)  
    
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
            
            # add the codons that were mutated when generating sgRNA 
            if(OverlapCodon1[1] != coords$codon[1]){
              PAM_mutated[[ID]][["sgRNA_mut_codon_overlap1"]] <- OverlapCodon1   
            }
            
            if(OverlapCodon2[1] != coords$codon[1]){
              PAM_mutated[[ID]][["sgRNA_mut_codon_overlap2"]] <- OverlapCodon2   
            }
            
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
    
    
    ############################################################
    # test if existing mutations generate restriction sites,
    # add them to the output list
    ############################################################
    
    # define the wild-type site assay for checking non-cutters
    wt_seq <- coords[["sequence"]]
    
    # subset the full sequence to that amplified by primers
    wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
    
    # load restriction enzyme data
    all_enzymes <- names(RE_SITES)
    
    # compute the non-cutter enzymes
    non_cutters <- getNonCutters(all_enzymes, wt_site_assay)    
    
    ######################################
    # iteration over the PAM_mutated list
    ######################################
    
    for(j in 1:length(PAM_mutated)){
      
      #################################################
      # extract relevant parts of the current item 
      #################################################
      
      # get full genomicString sequence and subset to the sequence spanned by primers
      mut_site_assay <- substr(PAM_mutated[[j]][["site_assay"]], forw_primer_pos[1], rev_primer_pos[2])
      
      # make a reverse complement of the mut_site_assay
      rc_mut_site_assay <- as.character(reverseComplement(DNAString(mut_site_assay)))
      
      # get coordinates of the codon in the current list item
      codon_coords <- PAM_mutated[[j]][["codon_coords"]] - offset
      codon_diffs_coords <- PAM_mutated[[j]][["codon_diffs_coords"]] - offset
      
      # make a flag variable to keep track of finding sites that 
      # have been introduced
      SITES_EXIST_FLAG <- FALSE
      
      # subset the site assay to a small 30-nt string for testing
      cur_test_string <- substr(mut_site_assay, codon_coords[1] - 15, codon_coords[1] + 18)
      
      # run a function to get all enzymes that cut 
      cutters <- getCutters(non_cutters, cur_test_string)
      
      # iterate over enzymes that cut the test string around the codon
      if(length(cutters) > 0){
        
        # initialize the list for enzyme sites
        PAM_mutated[[j]][["enzymes"]] <- list()
        enzID <- 1
        
        # add all of the cutter enzymes to a list inside an existing oligo entry
        for(cutterEnzyme in cutters){
          
          # get cut site
          cut_site <- as.character(RE_SITES[[cutterEnzyme]]$RE_site)
          
          # test if the cut site is in the forward strand
          if(str_detect(toupper(mut_site_assay), cut_site)){
            
            # get coordinates
            site_coords <- c(str_locate(toupper(mut_site_assay), cut_site))    
            
            # store site orientation 
            site_oriented <- "sense"
            
          }else{
            # get coordinates
            site_coords <- c(rev(nchar(rc_mut_site_assay) - str_locate(toupper(rc_mut_site_assay), cut_site) + 1))
            
            # store site orientation 
            site_oriented <- "anti"
          }
          
          # store the restriction site information
          PAM_mutated[[j]][["enzymes"]][[enzID]] <- list()
          
          PAM_mutated[[j]][["enzymes"]][[enzID]][["RE_enzyme"]] <- cutterEnzyme
          PAM_mutated[[j]][["enzymes"]][[enzID]][["RE_site"]] <- RE_SITES[[cutterEnzyme]]$Sequence
          
          # storing the restriction site coordinates in the original coordinates
          PAM_mutated[[j]][["enzymes"]][[enzID]][["RE_site_coords"]] <- site_coords + offset
          
          # store site orientation 
          PAM_mutated[[j]][["enzymes"]][[enzID]][["RE_site_oriented"]] <- site_oriented
          
          # update the enzID counter
          enzID <- enzID + 1
          
        } # end of cutters for loop
        
      } # end of if statement  for testing how many enzymes cut the sequence
      

    } 
    ##################################
    # end of PAM_mutated iteration
    ##################################
    
    PAM_mutated
  }) # end of PAM_mutations
  
  

  # adding silent mutations to introduce restriction sites
  # gets input from PAM_mutations()
  # this function will generate the final output to visualize the oligos and knock-in design
  
  REsite_silent_mutations <- reactive({
    
    # create a list for output
    output_REsites <- list()
    
    # main a counter variable
    ID <- 1
    
    ##############################################
    # coordinates extraction 
    ##############################################
    
    # get the coordinates of sgRNA and PAM which would be 
    # common to all mutations
    coords <-strategyCoords()
    
    codonPos <- coords[["codon"]]
    exonPos <- coords[["exon"]]
    exonStart <- exonPos[1]
    exonEnd <- exonPos[2]
    
    # sgRNA spacer
    sgRNA_pos <- coords[["sgRNA"]]
    
    # PAM
    pam_pos <- coords[["PAM"]]
    
    # extract the original sequence
    wt_seq <- coords[["sequence"]]
    
    # primers
    forw_primer_pos <- coords[["forw_primer"]]
    rev_primer_pos <- coords[["rev_primer"]] 
    
    ###########################################
    # conversion step
    ###########################################
    
    # subset the full sequence to that amplified by primers
    wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
    
    # offset all the positions
    offset <- forw_primer_pos[1] -1 # the length of the sequence before the primer
    
    codonPos <- codonPos - offset
    exonPos <- exonPos - offset
    exonStart <- exonStart - offset
    exonEnd <- exonEnd - offset
    
    # sgRNA spacer
    sgRNA_pos <- sgRNA_pos - offset
    
    # PAM
    pam_pos <- pam_pos - offset
    
    # primer positions
    forw_primer_pos <- forw_primer_pos - offset
    rev_primer_pos <- rev_primer_pos - offset
    
    ###############################################
    # find non-cutter enzymes in the wt_site_assay
    ###############################################
    
    # load restriction enzyme data
    all_enzymes <- names(RE_SITES)
    
    # compute the non-cutter enzymes
    non_cutters <- getNonCutters(all_enzymes, wt_site_assay)
    
    ##################################################################
    # main part of the function - processing previous designs 
    # to introduce restriction sites by synonymous mutations
    ##################################################################
    
    # get the previously mutated site assays to introduce either PAM or sgRNA mutations
    pam_muts <- PAM_mutations()
    
    for(item in pam_muts){
      
      #################################################
      # extract relevant parts of item and offset them
      # subset the sequence to the site assay
      # store coordinates of mutated codons
      #################################################
      
      # get full genomicString sequence and subset to the sequence spanned by primers
      # use the old coordinates since the sequence in item is the full genomicString
      # offset is added to get the old coordinates
      mut_site_assay <- substr(item[["site_assay"]], forw_primer_pos[1] + offset, rev_primer_pos[2] + offset)
      
      # make a reverse complement of the mut_site_assay
      rc_mut_site_assay <- as.character(reverseComplement(DNAString(mut_site_assay)))
      
      # copy the previous contents of the list from PAM_mutations() output
      new_codon <- item[["new_codon"]]
      AA_mut <- item[["AA_mutation"]]
      codon_coords <- item[["codon_coords"]] - offset
      codon_diffs_coords <- item[["codon_diffs_coords"]] - offset
      
      # start a list for all additional mutated codons
      mutated_codons <- list() # the list will contain the coordinate vectors of already mutated codons
      codon_id <- 1
      
      # check if PAM was mutated and add the coordinates to the list
      # store the other parameters of the mutation in variables
      if(item[["PAM_mutant_codon"]] != "none"){
        
        PAM_mutant_codon <- item[["PAM_mutant_codon"]]
        PAM_mut_codon_coords <- item[["PAM_mut_codon_coords"]] - offset # offset coordinates
        
        # store the mutated PAM codon in the mutated_codons list
        mutated_codons[[codon_id]] <- PAM_mut_codon_coords
        codon_id <- codon_id + 1
        
        # assign the sequence differences resulting from the PAM-related mutation
        # to a variable
        PAM_mut_codon_diffs <- item[["PAM_mut_codon_diffs"]] - offset
        
      } else{
        PAM_mutant_codon <- item[["PAM_mutant_codon"]]
      } # end of if-else statement for PAM mutation
      
      
      # check if there were any sgRNA mutations and store their parameters in variables
      if(item[["sgRNA_mutations"]]){
        
        # store all mutation difference
        sgRNA_mut_codon_diffs <- item[["sgRNA_mut_codon_diffs"]] - offset
        
        # check that and which overlap codons are present in the data  
        if( !is.null( item[["sgRNA_mut_codon_overlap1"]] ) ){
          sgRNA_mut_codon_overlap1 <- item[["sgRNA_mut_codon_overlap1"]] - offset
          
          mutated_codons[[codon_id]] <- sgRNA_mut_codon_overlap1
          codon_id <- codon_id + 1
        }
        
        
        if( !is.null( item[["sgRNA_mut_codon_overlap2"]] ) ){
          sgRNA_mut_codon_overlap2 <- item[["sgRNA_mut_codon_overlap2"]] - offset
          
          mutated_codons[[codon_id]] <- sgRNA_mut_codon_overlap2
          codon_id <- codon_id + 1
          
        }
        
      } # end of sgRNA mutations statements
      # end of mutate codon storage
      
      ############################################################
      # test if existing mutations generate restriction sites
      ############################################################
      
      # make a flag variable to keep track of finding sites that 
      # have been introduced
      SITES_EXIST_FLAG <- FALSE
      
      # subset the site assay to a small 30-nt string for testing
      cur_test_string <- substr(mut_site_assay, codon_coords[1] - 15, codon_coords[1] + 18)
      
      # run a function to get all enzymes that cut 
      cutters <- getCutters(non_cutters, cur_test_string)
      
      # iterate over enzymes that cut the test string around the codon
      if(length(cutters) > 0){
        
        # sites found, so update the flag variable
        SITES_EXIST_FLAG <- TRUE
        
        # iterate over each enzyme to get the site coordinates and store the results
        #########################################
        # STORAGE OF RESULTS TO THE OUTPUT LIST
        #########################################
        
        # initiate the list for this ID
        output_REsites[[ID]] <-  list()
        
        # store the current site assay
        output_REsites[[ID]][["site_assay"]] <- mut_site_assay
        
        # store data on the main codon mutation
        output_REsites[[ID]][["new_codon"]] <- new_codon
        output_REsites[[ID]][["AA_mutation"]] <- AA_mut
        output_REsites[[ID]][["codon_coords"]] <- codon_coords
        output_REsites[[ID]][["codon_diffs_coords"]] <- codon_diffs_coords
        
        # add the information on PAM mutations to the new output data structure
        if(item[["PAM_mutant_codon"]] != "none"){
          
          output_REsites[[ID]][["PAM_mutant_codon"]] <- PAM_mutant_codon
          output_REsites[[ID]][["PAM_mut_codon_coords"]] <- PAM_mut_codon_coords            
          output_REsites[[ID]][["PAM_mut_codon_diffs"]] <- PAM_mut_codon_diffs                  
          output_REsites[[ID]][["sgRNA_mutations"]] <- FALSE
          
        } else{
          
          output_REsites[[ID]][["PAM_mutant_codon"]] <- "none"
          
        } # end of if-else statement for PAM mutation
        
        # check if there were any sgRNA mutations and store their parameters in variables
        if(item[["sgRNA_mutations"]]){
          
          # store an indicator for sgRNA_mutations
          output_REsites[[ID]][["sgRNA_mutations"]] <- item[["sgRNA_mutations"]]
          
          # store all mutation difference
          output_REsites[[ID]][["sgRNA_mut_codon_diffs"]] <- sgRNA_mut_codon_diffs
          
          # check that and which overlap codons are present in the data  
          if( !is.null( item[["sgRNA_mut_codon_overlap1"]] ) ){
            output_REsites[[ID]][["sgRNA_mut_codon_overlap1"]] <- sgRNA_mut_codon_overlap1
          }
          
          if( !is.null( item[["sgRNA_mut_codon_overlap2"]] ) ){
            output_REsites[[ID]][["sgRNA_mut_codon_overlap2"]] <- sgRNA_mut_codon_overlap2
          }
          
        } else{
          
          # store an indicator for sgRNA_mutations
          output_REsites[[ID]][["sgRNA_mutations"]] <- item[["sgRNA_mutations"]]            
          
        } # end of sgRNA mutations statements                    
        
        
        # initialize the list for enzyme sites
        output_REsites[[ID]][["enzymes"]] <- list()
        enzID <- 1
        
        # add all of the cutter enzymes to a list inside an existing oligo entry
        for(cutterEnzyme in cutters){
          
          # get cut site
          cut_site <- as.character(RE_SITES[[cutterEnzyme]]$RE_site)
          
          # test if the cut site in the cut site is in the main ("forward") strand
          if(str_detect(toupper(mut_site_assay), cut_site)){
            
            # get coordinates
            site_coords <- c(str_locate(toupper(mut_site_assay), cut_site))    
            
            # store site orientation 
            site_oriented <- "sense"
            
          }else{
            # get coordinates
            site_coords <- c(rev(nchar(rc_mut_site_assay) - str_locate(toupper(rc_mut_site_assay), cut_site) + 1))
            
            # store site orientation 
            site_oriented <- "anti"
            
          }
          
          # store the restriction site information
          output_REsites[[ID]][["enzymes"]][[enzID]] <- list()
          
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_enzyme"]] <- cutterEnzyme
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site"]] <- RE_SITES[[cutterEnzyme]]$Sequence
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_coords"]] <- site_coords
          
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_oriented"]] <- site_oriented
          
          # update the enzID counter
          enzID <- enzID + 1
          
        } # end of cutters for loop
        
        # update the main ID variable
        ID <- ID + 1
        
      } # end of if statement  for testing how many enzymes cut the sequence
      

      ###############################################
      # existing mutations do not introduce any sites

        # flag for introducing a site by a synonymous mutation
        REsite_BY_MUTATION <- FALSE
        
        #############################################################
        # codon selection for mutating to introduce restriction sites
        #############################################################
        
        ################################################################
        # The idea: test up to 3 codons on each side of the target codon
        # taking each time a codon on the left side or right side
        # test them for NOT being inside the mutated codon list and for being WITHIN the exon,
        # update the respective codon pointers and if the test was positive, add the codons to the list
        # for the candidates to be mutated
        
        
        # list for codons to be mutated
        codons2mut4REsites <- list()
        sel_codon_id <- 1
        
        # total number 
        num_checked <- 0
        
        # make pointers
        left_pointer <- 1
        right_pointer <- 1
        
        
        # loop to iterate over potential codons
        while(length(codons2mut4REsites) <= 3 & num_checked < 6 ){
          
          # select a codon 5' (left) from the target codon
          leftCodon <- codonPos - left_pointer*3
          
          # update num_checked variable
          num_checked <- num_checked + 1
          left_pointer <- left_pointer + 1
          
          # check if this codon is OK
          if( !vector_in_list(mutated_codons, leftCodon) & (leftCodon[1] >=  exonStart) ){
            
            codons2mut4REsites[[sel_codon_id]] <- leftCodon
            sel_codon_id <- sel_codon_id + 1
            
          }          
          
          
          # select a codon 3' (right) from the target codon
          rightCodon <- codonPos + right_pointer*3
          
          # update num_checked variable and the pointer
          num_checked <- num_checked + 1
          right_pointer <- right_pointer + 1
          
          # check if this codon is OK
          if( !vector_in_list(mutated_codons, rightCodon) & (rightCodon[1] <=  exonEnd) ){
            
            codons2mut4REsites[[sel_codon_id]] <- rightCodon
            sel_codon_id <- sel_codon_id + 1
            
          }          
          
          
        } # end of while loop for finding codons to be mutated
        
        
        ###########################################################
        # synonymous mutations in the selected codons
        # and evaluation of non-cutters on the new mutant sequences
        ###########################################################
        
        # define a new vector of non-cutter enzymes
        # the idea is that we want to remove all enzymes that have already been found
        # and focus on the newly found ones
        non_cutters_mut <- getNonCutters(all_enzymes, mut_site_assay)
        
        
        # iterate over all potential codons to be selected
        for(codon in codons2mut4REsites){
          
          # find the non-identical synonymous codons
          selCodonSeq = toupper(substr(mut_site_assay, codon[1], codon[2]))
          
          # get all possible codons for the encoded amino acid
          aa_codons = REV_GENETIC_CODE[[ GENETIC_CODE[[selCodonSeq]] ]]
          
          # get codons that are not identical to the current codon
          codons_not_same = aa_codons[aa_codons != selCodonSeq]
          
          # perform all possible synonymous codon replacements
          for(codon_nonID in codons_not_same){
            
            # make a mutant site assay version, store in a temp variable not to interfere with subsequent steps
            mut_cds <- paste(substr(mut_site_assay, 1, codon[1]-1), codon_nonID, substr(mut_site_assay, codon[2] + 1, nchar(mut_site_assay)), sep = "")
            
            # make a reverse complement of the mut_cds
            rc_mut_cds <- as.character(reverseComplement(DNAString(mut_cds)))
            
            # subset the site assay to a small 30-nt string for testing
            cur_test_string <- substr(mut_cds, codon[1] - 10, codon[1] + 13)
            
            # run a function to get all enzymes that cut 
            new_cutters <- getCutters(non_cutters_mut, cur_test_string)
            
            #####################################
            # CUT IS SUCCESSFUL
            #####################################
            
            if(length(new_cutters) > 0){
              
              # update the flag variable
              REsite_BY_MUTATION <- TRUE
              
              #########################################
              # STORAGE OF RESULTS TO THE OUTPUT LIST
              #########################################
              
              # initiate the list for this ID
              output_REsites[[ID]] = list()
              
              # store the current site assay
              output_REsites[[ID]][["site_assay"]] <- mut_cds
              
              # store data on the main codon mutation
              output_REsites[[ID]][["new_codon"]] <- new_codon
              output_REsites[[ID]][["AA_mutation"]] <- AA_mut
              output_REsites[[ID]][["codon_coords"]] <- codon_coords
              output_REsites[[ID]][["codon_diffs_coords"]] <- codon_diffs_coords
              
              # add the information on PAM mutations to the new output data structure
              if(item[["PAM_mutant_codon"]] != "none"){
                
                output_REsites[[ID]][["PAM_mutant_codon"]] <- PAM_mutant_codon
                output_REsites[[ID]][["PAM_mut_codon_coords"]] <- PAM_mut_codon_coords            
                output_REsites[[ID]][["PAM_mut_codon_diffs"]] <- PAM_mut_codon_diffs                  
                output_REsites[[ID]][["sgRNA_mutations"]] <- FALSE
                
              } else{
                
                output_REsites[[ID]][["PAM_mutant_codon"]] <- "none"
                
              } # end of if-else statement for PAM mutation
              
              # check if there were any sgRNA mutations and store their parameters in variables
              if(item[["sgRNA_mutations"]]){
                
                # store the sgRNA_mutations indicator
                output_REsites[[ID]][["sgRNA_mutations"]] <- item[["sgRNA_mutations"]]
                
                # store all mutation difference
                output_REsites[[ID]][["sgRNA_mut_codon_diffs"]] <- sgRNA_mut_codon_diffs
                
                # check that and which overlap codons are present in the data  
                if( !is.null( item[["sgRNA_mut_codon_overlap1"]] ) ){
                  output_REsites[[ID]][["sgRNA_mut_codon_overlap1"]] <- sgRNA_mut_codon_overlap1
                }
                
                if( !is.null( item[["sgRNA_mut_codon_overlap2"]] ) ){
                  output_REsites[[ID]][["sgRNA_mut_codon_overlap2"]] <- sgRNA_mut_codon_overlap2
                }
                
              } else{
                # store the sgRNA_mutations indicator
                output_REsites[[ID]][["sgRNA_mutations"]] <- item[["sgRNA_mutations"]]
                
              } # end of sgRNA mutations statements                    
              
              # initialize the list for enzyme sites
              output_REsites[[ID]][["enzymes"]] <- list()
              enzID <- 1
              
              # define a new full list of cutter enzymes inclusing new ones by synonymous 
              # mutations and the original PAM/sgRNA ones
              all_cutters <- c(new_cutters, cutters)
              
              
              # iterate over each enzyme to get the site coordinates and store the results
              for(cutterEnzyme in all_cutters){
                
                # get cut site
                cut_site <- as.character(RE_SITES[[cutterEnzyme]]$RE_site)
                
                # initialize an empty vector for coordinates
                site_coords <- c()
                
                # test if the cut site in the cut site is in the main ("forward") strand
                if(str_detect(toupper(mut_cds), cut_site)){
                  
                  # get coordinates
                  site_coords <- c(str_locate(toupper(mut_cds), cut_site))
                  
                  # store site orientation
                  site_oriented <- "sense"
                  
                }else if(str_detect(toupper(rc_mut_cds), cut_site)){
                  # get coordinates
                  site_coords <- c(rev(nchar(rc_mut_cds) - str_locate(toupper(rc_mut_cds), cut_site) + 1 ))
                  
                  # store site orientation
                  site_oriented <- "anti"
                } 
                
                
                # condition the output on the correct matching of a restriction site
                if(length(site_coords) > 0){

                  # store the restriction site information
                  output_REsites[[ID]][["enzymes"]][[enzID]] <- list()
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_enzyme"]] <- cutterEnzyme
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site"]] <- RE_SITES[[cutterEnzyme]]$Sequence
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_coords"]] <- site_coords
                  
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_oriented"]] <- site_oriented
                  
                  # store sequence changes that led to the site introduction
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_codon"]] <- codon_nonID 
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_codon_coords"]] <- codon
                  output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_codon_diffs"]] <- getCoordinatesDiffs(codon_nonID, selCodonSeq, codon[1])   
                  
                  enzID <- enzID + 1
                  
                }
                
              } # end of cutters for loop
              
              # update ID
              ID <- ID + 1
              
            } # end of if statement  for testing how many enzymes cut the sequence
            
          } # end of the for loop over synonymous codons
          
        } # end of for loop over codons2mut4REsites
        
      # } # end of if block the case SITES_EXIST_FLAG is FALSE
      
      ###############################################################################
      # Evaluation if the current item was mutated to introduce any restriction sites
      ###############################################################################
      
      # the condition below is TRUE when no sites have been introduced by mutations
      if(!REsite_BY_MUTATION ){
        
        #########################################
        # STORAGE OF RESULTS TO THE OUTPUT LIST
        #########################################
        
        # initiate the list for this ID
        output_REsites[[ID]] = list()
        
        # store the current site assay
        output_REsites[[ID]][["site_assay"]] <- mut_site_assay
        
        # store data on the main codon mutation
        output_REsites[[ID]][["new_codon"]] <- new_codon
        output_REsites[[ID]][["AA_mutation"]] <- AA_mut
        output_REsites[[ID]][["codon_coords"]] <- codon_coords
        output_REsites[[ID]][["codon_diffs_coords"]] <- codon_diffs_coords
        
        # add the information on PAM mutations to the new output data structure
        if(item[["PAM_mutant_codon"]] != "none"){
          
          output_REsites[[ID]][["PAM_mutant_codon"]] <- PAM_mutant_codon
          output_REsites[[ID]][["PAM_mut_codon_coords"]] <- PAM_mut_codon_coords            
          output_REsites[[ID]][["PAM_mut_codon_diffs"]] <- PAM_mut_codon_diffs                  
          output_REsites[[ID]][["sgRNA_mutations"]] <- FALSE
          
        } else{
          
          output_REsites[[ID]][["PAM_mutant_codon"]] <- "none"
          
        } # end of if-else statement for PAM mutation
        
        # check if there were any sgRNA mutations and store their parameters in variables
        if(item[["sgRNA_mutations"]]){
          
          # store an indicator for sgRNA_mutations
          output_REsites[[ID]][["sgRNA_mutations"]] <- item[["sgRNA_mutations"]]
          
          # store all mutation difference
          output_REsites[[ID]][["sgRNA_mut_codon_diffs"]] <- sgRNA_mut_codon_diffs
          
          # check that and which overlap codons are present in the data  
          if( !is.null( item[["sgRNA_mut_codon_overlap1"]] ) ){
            output_REsites[[ID]][["sgRNA_mut_codon_overlap1"]] <- sgRNA_mut_codon_overlap1
          }
          
          if( !is.null( item[["sgRNA_mut_codon_overlap2"]] ) ){
            output_REsites[[ID]][["sgRNA_mut_codon_overlap2"]] <- sgRNA_mut_codon_overlap2
          }
          
        } else{
          
          # store an indicator for sgRNA_mutations
          output_REsites[[ID]][["sgRNA_mutations"]] <- item[["sgRNA_mutations"]]            
          
        } # end of sgRNA mutations statements                    
        
        # indicate that no restriction sites have been found
        output_REsites[[ID]][["RE_site"]] <- "none"
        
        # update ID
        ID <- ID + 1
        
      } # end of if statement for REsite_BY_MUTATION flag variable test 
      
    }# end of for loop over PAM_mutations designs
    
    # return the output data structure
    return(output_REsites)
  }) # end of REsite_silent_mutations
  
  
  noPAMmuts_REsiteSilentMuts <- reactive({
    
    # create a list for output
    output_REsites <- list()
    
    # main a counter variable
    ID <- 1
    
    ##############################################
    # coordinates extraction 
    ##############################################
    
    # get the coordinates of sgRNA and PAM which would be 
    # common to all mutations
    coords <-strategyCoords()
    
    codonPos <- coords[["codon"]]
    exonPos <- coords[["exon"]]
    exonStart <- exonPos[1]
    exonEnd <- exonPos[2]
    
    # sgRNA spacer
    sgRNA_pos <- coords[["sgRNA"]]
    
    # PAM
    pam_pos <- coords[["PAM"]]
    
    # extract the original sequence
    wt_seq <- coords[["sequence"]]
    
    # primers
    forw_primer_pos <- coords[["forw_primer"]]
    rev_primer_pos <- coords[["rev_primer"]] 
    
    ###########################################
    # conversion step
    ###########################################
    
    # subset the full sequence to that amplified by primers
    wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
    
    # offset all the positions
    offset <- forw_primer_pos[1] -1 # the length of the sequence before the primer
    
    codonPos <- codonPos - offset
    exonPos <- exonPos - offset
    exonStart <- exonStart - offset
    exonEnd <- exonEnd - offset
    
    # sgRNA spacer
    sgRNA_pos <- sgRNA_pos - offset
    
    # PAM
    pam_pos <- pam_pos - offset
    
    # primer positions
    forw_primer_pos <- forw_primer_pos - offset
    rev_primer_pos <- rev_primer_pos - offset
    
    ###############################################
    # find non-cutter enzymes in the wt_site_assay
    ###############################################
    
    # load restriction enzyme data
    all_enzymes <- names(RE_SITES)
    
    # compute the non-cutter enzymes
    non_cutters <- getNonCutters(all_enzymes, wt_site_assay)
    
    
    ##################################################################
    # main part of the function - processing previous designs 
    # to introduce restriction sites by synonymous mutations
    ##################################################################
    
    # get the previously mutated site assays where codon mutations were introduced 
    codon_muts <- codonMutations()
    
    for(item in codon_muts){
      
      #################################################
      # extract relevant parts of item and offset them
      # subset the sequence to the site assay
      # store coordinates of mutated codons
      #################################################
      
      # get full genomicString sequence and subset to the sequence spanned by primers
      # use the old coordinates since the sequence in item is the full genomicString
      # offset is added to get the old coordinates
      mut_site_assay <- substr(item[["site_assay"]], forw_primer_pos[1] + offset, rev_primer_pos[2] + offset)
      
      # make a reverse complement of the mut_site_assay
      rc_mut_site_assay <- as.character(reverseComplement(DNAString(mut_site_assay)))
      
      # copy the previous contents of the list from codonMutations() output
      new_codon <- item[["new_codon"]]
      AA_mut <- item[["AA_mutation"]]
      codon_coords <- item[["codon_coords"]] - offset
      codon_diffs_coords <- item[["codon_diffs_coords"]] - offset
      
      # start a list for all additional mutated codons
      mutated_codons <- list() # the list will contain the coordinate vectors of already mutated codons
      codon_id <- 1
      
      
      ############################################################
      # test if existing mutations generate restriction sites
      ############################################################
      
      # make a flag variable to keep track of finding sites that 
      # have been introduced
      SITES_EXIST_FLAG <- FALSE
      
      # subset the site assay to a small 30-nt string for testing
      cur_test_string <- substr(mut_site_assay, codon_coords[1] - 15, codon_coords[1] + 18)
      
      # run a function to get all enzymes that cut 
      cutters <- getCutters(non_cutters, cur_test_string)
      
      # iterate over enzymes that cut the test string around the codon
      if(length(cutters) > 0){
        
        # sites found, so update the flag variable
        SITES_EXIST_FLAG <- TRUE
        
        # iterate over each enzyme to get the site coordinates and store the results
        #########################################
        # STORAGE OF RESULTS TO THE OUTPUT LIST
        #########################################
        
        # initiate the list for this ID
        output_REsites[[ID]] <-  list()
        
        # store the current site assay
        output_REsites[[ID]][["site_assay"]] <- mut_site_assay
        
        # store data on the main codon mutation
        output_REsites[[ID]][["new_codon"]] <- new_codon
        output_REsites[[ID]][["AA_mutation"]] <- AA_mut
        output_REsites[[ID]][["codon_coords"]] <- codon_coords
        output_REsites[[ID]][["codon_diffs_coords"]] <- codon_diffs_coords


        # initialize the list for enzyme sites
        output_REsites[[ID]][["enzymes"]] <- list()
        enzID <- 1
        
        # add all of the cutter enzymes to a list inside an existing oligo entry
        for(cutterEnzyme in cutters){
          
          # get cut site
          cut_site <- as.character(RE_SITES[[cutterEnzyme]]$RE_site)
          
          # test if the cut site in the cut site is in the main ("forward") strand
          if(str_detect(toupper(mut_site_assay), cut_site)){
            
            # get coordinates
            site_coords <- c(str_locate(toupper(mut_site_assay), cut_site))    
    
            # store site orientation 
            site_oriented <- "sense"
            
          }else{
            # get coordinates
            site_coords <- c(rev(nchar(rc_mut_site_assay) - str_locate(toupper(rc_mut_site_assay), cut_site) + 1))
            
            # store site orientation 
            site_oriented <- "anti"
          }
          
          # store the restriction site information
          output_REsites[[ID]][["enzymes"]][[enzID]] <- list()
          
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_enzyme"]] <- cutterEnzyme
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site"]] <- RE_SITES[[cutterEnzyme]]$Sequence
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_coords"]] <- site_coords
          
          output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_oriented"]] <- site_oriented
          
          # update the enzID counter
          enzID <- enzID + 1
          
        } # end of cutters for loop
        
        # update the main ID variable
        ID <- ID + 1
        
      } # end of if statement  for testing how many enzymes cut the sequence
      
      ###############################################
      # existing mutations do not introduce any sites

      
        # flag for introducing a site by a synonymous mutation
        REsite_BY_MUTATION <- FALSE
        
        #############################################################
        # codon selection for mutating to introduce restriction sites
        #############################################################
        
        ################################################################
        # The idea: test up to 3 codons on each side of the target codon
        # taking each time a codon on the left side or right side
        # test them for NOT being inside the mutated codon list and for being WITHIN the exon,
        # update the respective codon pointers and if the test was positive, add the codons to the list
        # for the candidates to be mutated
        
        
        # list for codons to be mutated
        codons2mut4REsites <- list()
        sel_codon_id <- 1
        
        # total number 
        num_checked <- 0
        
        # make pointers
        left_pointer <- 1
        right_pointer <- 1
        
        
        # loop to iterate over potential codons
        while(length(codons2mut4REsites) <= 3 & num_checked < 6 ){
          
          # select a codon 5' (left) from the target codon
          leftCodon <- codonPos - left_pointer*3
          
          # update num_checked variable
          num_checked <- num_checked + 1
          left_pointer <- left_pointer + 1
          
          # check if this codon is OK
          if( !vector_in_list(mutated_codons, leftCodon) & (leftCodon[1] >=  exonStart) ){
            
            codons2mut4REsites[[sel_codon_id]] <- leftCodon
            sel_codon_id <- sel_codon_id + 1
            
          }          
          
          
          # select a codon 3' (right) from the target codon
          rightCodon <- codonPos + right_pointer*3
          
          # update num_checked variable and the pointer
          num_checked <- num_checked + 1
          right_pointer <- right_pointer + 1
          
          # check if this codon is OK
          if( !vector_in_list(mutated_codons, rightCodon) & (rightCodon[1] <=  exonEnd) ){
            
            codons2mut4REsites[[sel_codon_id]] <- rightCodon
            sel_codon_id <- sel_codon_id + 1
            
          }          
          
          
        } # end of while loop for finding codons to be mutated
        
        
        ###########################################################
        # synonymous mutations in the selected codons
        # and evaluation of non-cutters on the new mutant sequences
        ###########################################################
        
        # define a new vector of non-cutter enzymes
        # the idea is that we want to remove all enzymes that have already been found
        # and focus on the newly found ones
        non_cutters_mut <- getNonCutters(all_enzymes, mut_site_assay)
        
        # iterate over all potential codons to be selected
        for(codon in codons2mut4REsites){
          
          # find the non-identical synonymous codons
          selCodonSeq = toupper(substr(mut_site_assay, codon[1], codon[2]))
          
          # get all possible codons for the encoded amino acid
          aa_codons = REV_GENETIC_CODE[[ GENETIC_CODE[[selCodonSeq]] ]]
          
          # get codons that are not identical to the current codon
          codons_not_same = aa_codons[aa_codons != selCodonSeq]
          
          # perform all possible synonymous codon replacements
          for(codon_nonID in codons_not_same){
            
            # make a mutant site assay version, store in a temp variable not to interfere with subsequent steps
            mut_cds <- paste(substr(mut_site_assay, 1, codon[1]-1), codon_nonID, substr(mut_site_assay, codon[2] + 1, nchar(mut_site_assay)), sep = "")
            
            # make a reverse complement of the mut_cds
            rc_mut_cds <- as.character(reverseComplement(DNAString(mut_cds)))
            
            # subset the site assay to a small 30-nt string for testing
            cur_test_string <- substr(mut_cds, codon[1] - 10, codon[1] + 13)
            
            # run a function to get all enzymes that cut 
            new_cutters <- getCutters(non_cutters_mut, cur_test_string)
            
            #####################################
            # CUT IS SUCCESSFUL
            #####################################
            
            if(length(new_cutters) > 0){
              
              # update the flag variable
              REsite_BY_MUTATION <- TRUE
              
              #########################################
              # STORAGE OF RESULTS TO THE OUTPUT LIST
              #########################################
              
              # initiate the list for this ID
              output_REsites[[ID]] = list()
              
              # store the current site assay
              output_REsites[[ID]][["site_assay"]] <- mut_cds
              
              # store data on the main codon mutation
              output_REsites[[ID]][["new_codon"]] <- new_codon
              output_REsites[[ID]][["AA_mutation"]] <- AA_mut
              output_REsites[[ID]][["codon_coords"]] <- codon_coords
              output_REsites[[ID]][["codon_diffs_coords"]] <- codon_diffs_coords
              

              # initialize the list for enzyme sites
              output_REsites[[ID]][["enzymes"]] <- list()
              enzID <- 1
              
              # define a new full list of cutter enzymes inclusing new ones by synonymous 
              # mutations and the original PAM/sgRNA ones
              all_cutters <- c(new_cutters, cutters)
              
              # iterate over each enzyme to get the site coordinates and store the results
              for(cutterEnzyme in all_cutters){
                
                # get cut site
                cut_site <- as.character(RE_SITES[[cutterEnzyme]]$RE_site)
                
                # test if the cut site in the cut site is in the main ("forward") strand
                if(str_detect(toupper(mut_cds), cut_site)){
                  
                  # get coordinates
                  site_coords <- c(str_locate(toupper(mut_cds), cut_site))
                  
                  # store site orientation
                  site_oriented <- "sense"
                  
                }else if(str_detect(toupper(rc_mut_cds), cut_site)){
                  # get coordinates
                  site_coords <- c(rev(nchar(rc_mut_cds) - str_locate(toupper(rc_mut_cds), cut_site) + 1 ))
                  
                  # store site orientation
                  site_oriented <- "anti"
                  
                }
                
                if(length(site_coords) > 0){
                
                # store the restriction site information
                output_REsites[[ID]][["enzymes"]][[enzID]] <- list()
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_enzyme"]] <- cutterEnzyme
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site"]] <- RE_SITES[[cutterEnzyme]]$Sequence
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_coords"]] <- site_coords
                
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_oriented"]] <- site_oriented
                
                # store sequence changes that led to the site introduction
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_codon"]] <- codon_nonID 
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_codon_coords"]] <- codon
                output_REsites[[ID]][["enzymes"]][[enzID]][["RE_site_codon_diffs"]] <- getCoordinatesDiffs(codon_nonID, selCodonSeq, codon[1])   
                
                enzID <- enzID + 1
                
                } # end of if-statement for site coordinates
                
              } # end of cutters for loop
              
              # update ID
              ID <- ID + 1
              
            } # end of if statement  for testing how many enzymes cut the sequence
            
          } # end of the for loop over synonymous codons
          
        } # end of for loop over codons2mut4REsites
        
       # } # end of if block the case SITES_EXIST_FLAG is FALSE
      
      ###############################################################################
      # Evaluation if the current item was mutated to introduce any restriction sites
      ###############################################################################
      
      # the condition below is TRUE when no sites have been introduced by mutations
      if(!REsite_BY_MUTATION ){
        
        #########################################
        # STORAGE OF RESULTS TO THE OUTPUT LIST
        #########################################
        
        # initiate the list for this ID
        output_REsites[[ID]] = list()
        
        # store the current site assay
        output_REsites[[ID]][["site_assay"]] <- mut_site_assay
        
        # store data on the main codon mutation
        output_REsites[[ID]][["new_codon"]] <- new_codon
        output_REsites[[ID]][["AA_mutation"]] <- AA_mut
        output_REsites[[ID]][["codon_coords"]] <- codon_coords
        output_REsites[[ID]][["codon_diffs_coords"]] <- codon_diffs_coords
        
        # indicate that no restriction sites have been found
        output_REsites[[ID]][["RE_site"]] <- "none"
        
        # update ID
        ID <- ID + 1
        
      } # end of if statement for REsite_BY_MUTATION flag variable test 
      
    }# end of for loop over codonMutations designs
    
    # return the output data structure
    return(output_REsites)
  }) # end of noPAMmuts_REsiteSilentMuts
  
  ########################################################################
  # FINAL OUTPUT STAGE
  ########################################################################
  
  # The UI function will have to run the knockinDesign function to produce the complete description for all 
  # the designs that are possible with the current input 
  # it will then take the output of the knockinDesign and present all the individual designs
  
  # main handler to run oligo design
  observeEvent(input$run, {
    
    shinyjs::runjs("window.scrollTo(0,0);") 
    
    
    # a function to output oligos with introduced restriction sites
    output$finalOligos <- renderUI({
    
      # obtain the data on the overall strategy
      # make sure the reactive code only runs when you press "Submit" button
      isolate({ coords <- strategyCoords() })
      
      ###################################
      # outputStrategy function
      ###################################
      
      outputStrategy <- function(coords){
        
        # initialize the string for HTML output
        outputHTML <- '<br/><button type="button" class="btn btn-primary" style="width: 100%; font-size: 20px">Targeting strategy outline:</button><div class="jumbotron", style="width: 100%; word-wrap:break-word; display:inline-block; font-size: 16px;">'
        
        # process the input data to add codes to the codon, sgRNA and PAM positions
        # the primer and intervening positions can be added in a regular way
        
        for_primer_pos <- coords[["forw_primer"]]
        rev_primer_pos <- coords[["rev_primer"]]
        
        sequence = coords[["sequence"]]
        
        
        # generate the vectors of all important positions
        codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
        pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
        sgRNA_pos <- coords[["sgRNA"]][1]:coords[["sgRNA"]][2]
        exon_pos <- coords[["exon"]][1]: coords[["exon"]][2]
        all_special_pos <- sort(unique(c(codon_pos, pam_pos, sgRNA_pos)))
        
        # an improved algorithm:
        # 1. iterate from the start of forward primer to the end of reverse primer
        # 2. make EXON sequence bold and add relevant formatting to other special parts
        # of the sequence, keep intron letters regular lowercase.
        # 3. Remove the other loops or bulk substring commands
        
        
        # amplicon iteration
        for(i in for_primer_pos[1]:rev_primer_pos[2]){
          
          # forward primer
          if(i >= for_primer_pos[1] & i <= for_primer_pos[2]){
            
            # format according to whether a nucleotide is inside an exon
            if(i %in% exon_pos){
              outputHTML <- paste(outputHTML, "<strong><font style='BACKGROUND-COLOR: #90ee90; color: black'>",
                                  substr(sequence, i, i),"</font></strong>",sep = "" )              
            }else{
              outputHTML <- paste(outputHTML, "<font style='BACKGROUND-COLOR: #90ee90; color: black'>", substr(sequence, i, i),"</font>",sep = "" )              
              
            }
            next
          } # end forward primer
          
          
          # rev primer
          if(i >= rev_primer_pos[1] & i <= rev_primer_pos[2]){
            
            # format according to whether a nucleotide is inside an exon
            if(i %in% exon_pos){
              outputHTML <- paste(outputHTML, "<strong><font style='BACKGROUND-COLOR: #90ee90; color: black'>",
                                  substr(sequence, i, i),"</font></strong>",sep = "" )              
            }else{
              outputHTML <- paste(outputHTML, "<font style='BACKGROUND-COLOR: #90ee90; color: black'>", substr(sequence, i, i),"</font>",sep = "" )              
              
            }
            next
          }
          
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
              
              # make bold if the position is inside a codon
              if(i %in% exon_pos){
                # underlined + blue text
                outputHTML <- paste(outputHTML,"<u><strong><font style='color: #000080'>",
                                    substr(sequence, i,i), "</font></strong></u>", sep="")                
              }else{
                # underlined + blue text
                outputHTML <- paste(outputHTML,"<u><font style='color: #000080'>",
                                    substr(sequence, i,i), "</font></u>", sep="")
              }
              
            } #end of PAM NOT codon
            
            # sgRNA NOT codon
            if( (i %in% sgRNA_pos) && !(i %in% codon_pos) && !(i %in% pam_pos)){
              
              # make bold if the position is inside a codon
              if(i %in% exon_pos){
                # fuchsia-colored font
                outputHTML <- paste(outputHTML,"<strong><font style='color: fuchsia'>",
                                    substr(sequence, i,i), "</font></strong>", sep="")
              }else{
                # fuchsia-colored font
                outputHTML <- paste(outputHTML,"<font style='color: fuchsia'>",
                                    substr(sequence, i,i), "</font>", sep="")
              }
            } # end - sgRNA NOT codon
            
            
          } else{ # position is not labeled
            
            # make bold if the position is inside a codon
            if(i %in% exon_pos){
              outputHTML <- paste(outputHTML,"<strong>", substr(sequence, i,i), "</strong>", sep="")
            }else{
              outputHTML <- paste(outputHTML,substr(sequence, i,i), sep="")              
            }
            
          } # end of the if statements
          
        } # end for loop - amplicon iteration
        
        # add legend and finish this container
        outputHTML <- paste(outputHTML, "<br/>","<br/>", "<p style='font-size: 12pt; font-weight:bold'>Legend: ", 
                            "<strong><font style='BACKGROUND-COLOR: #90ee90; color: black'>", 
                            "primer","</font></strong>","  ", 
                            "<strong><font style='BACKGROUND-COLOR: yellow'>", "Codon",
                            "</font></strong>","  ",
                            "<u><strong><font style='color: #000080'>", "PAM sequence",
                            "</font></strong></u>","  ",
                            "<strong><font style='color: fuchsia'>", "sgRNA sequence",
                            "</font></strong>","  ", 
                            "<font style='font-size: 12pt;'>EXON  ","</font>", 
                            "<font style='font-size: 12pt; font-weight: normal'>intron", "</font>",
                            "</p>",
                            "</div>", sep = "")
        
        # add the header for the next output
        outputHTML <- paste(outputHTML, '<button type="button" class="btn btn-primary" style="width: 100%; font-size: 20px">Results of the oligo design:</button>',  sep = "") 
        
        outputHTML
      }
      
      
      ###################################
      # END outputStrategy function
      ###################################
      
      # generate the vectors of all important positions
      codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
      pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
      
      # define the original codon
      orig_codon <- substr(coords[["sequence"]], coords[["codon"]][1], coords[["codon"]][2] )
      
      # offset the coordinates for PAM and codon
      forw_primer_pos <- coords[["forw_primer"]]
      offset <- forw_primer_pos[1] - 1
      
      codon_pos <- codon_pos - offset
      pam_pos <- pam_pos - offset
      
      
      # initialize the output HTML
      strategyHTML <- outputStrategy(coords)
      
      ################################################
      # OUTPUT based on options: mutate_PAM and REsites
      ################################################
      
      # "NGG","NRG", "NGA", "NGCG", "NNGRRT", "NGG"
      # "TTTV", "TTTN", "TTN", "YTN", "NTTN", "YTTN"
      
  if(input$PAM %in% c("NGG","NRG", "NGA", "NGCG", "NNGRRT", "NGG")){
    
    # CASE 1: mutatePAM == "yes" AND "REsites" == "yes"
    if ( input$mutatePAM == "yes" & input$REsites == "yes" ) {
      
      # the output for complete mutation design including both PAM/sgRNA and RE site mutations
      isolate({ outputList <- REsite_silent_mutations() })
      
      # iterate over the list items
      lapply(1:length(outputList), function(j) {
        
        # get sequence
        sequence <- toupper(outputList[[j]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        PAM_muts <- c()
        sgRNA_muts <- c()
        REsite_muts <- c()
      
        # get new replacement codon sequence
        new_codon <- substr(sequence, codon_pos[1], codon_pos[3])
        
        # codon mutations
        codon_muts <- outputList[[j]][["codon_diffs_coords"]]
        
        # get PAM mutations is they are available
        if(outputList[[j]][["PAM_mutant_codon"]] != "none"){
          
          PAM_muts <- outputList[[j]][["PAM_mut_codon_diffs"]]
          
        }
        
        # get sgRNA mutations
        if(outputList[[j]][["sgRNA_mutations"]]){
          
          sgRNA_muts <- outputList[[j]][["sgRNA_mut_codon_diffs"]] 
          
        }     
        
        # get the codon mutation that introduced a restriction site
        if( "enzymes" %in% names(outputList[[j]])){
          
          # test if "RE_site_codon_diffs" in the relevant list
          if("RE_site_codon_diffs" %in% names(outputList[[j]][["enzymes"]][[1]])){
            REsite_muts <- outputList[[j]][["enzymes"]][[1]][["RE_site_codon_diffs"]]
          }  
        } 
        
        
        # combine all mutations
        mutated <- c(codon_muts, PAM_muts, sgRNA_muts, REsite_muts)
        
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
        
        
        # define the start and end of oligos for the purposes of correct output
        if(input$oriented == "sense"){
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[1] - 3 - input$leftArmLength
          oligoEnd <- pam_pos[1] - 3 + input$rightArmLength - 1
          
        }else{
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[3] + 3 - input$leftArmLength + 1
          oligoEnd <- pam_pos[3] + 3 + input$rightArmLength 
          
        }
        
        # organize all mutated positions into a single vector
        all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        # consider checking whether start and end of the oligo are less and more than any labeled positions
        # in the sequence and then updating them accordingly
        # Also check that neither of the coordinates is negative or beyond the sequence length,
        # change them accordingly
        
        # add the strategy HTML if the program is at the beginning of the lapply loop
        if(j == 1){
          # initialize the output HTML
          outputHTML <- HTML(paste(strategyHTML, "<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                              oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        } else{
          
          # initialize the output HTML
          outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        }
        ############################################ formatOligo function ##############################################
        
        outputHTML <- paste(outputHTML,  formatOligo(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated), sep="")
        
        ###########################################################################################
        
        # get all information on the relevant restriction enzyme site
        if("enzymes" %in% names(outputList[[j]])){
          
          # iterate over all known enzymes for that sequence
          for(enzymeData in outputList[[j]][["enzymes"]]){
            
            # store the restriction site information
            enzyme <- enzymeData[["RE_enzyme"]]
            site <- enzymeData[["RE_site"]]
            site_coords <- enzymeData[["RE_site_coords"]]
            site_oriented <- enzymeData[["RE_site_oriented"]]
            
            # change the coordinates of the restriction site according to the oligo orientation
            if(input$orientedOligo == "anti"){
              site_coords <- nchar(sequence) - rev(site_coords) + 1
            }
            
            
            # add a second line for the restriction site
            outputHTML <- paste(outputHTML, "<span style='opacity: 0;'>", substr(sequence, oligoStart, site_coords[1]-1), "</span>", 
                                "<strong>",substr(sequence, site_coords[1], site_coords[2]), "</strong>", " - ", "<span style='color:blue'>",
                                enzyme, "</span>", " (", "<span style=' font-weight: bold'>",  site,  "</span>", ") ",
                                calculateFragments(RE_SITES, enzyme, enzymeData[["RE_site_coords"]], site_oriented, sequence),"<br/>", sep = "")            
          }
          
          
        }
        
        ############################################ AS-PCR primers output #############################################
        
        # steps for outputting tables of primers and their characteristics
        
        ## 1. collect variables for AS-PCR primer designs
        forw_primer_pos <- coords[["forw_primer"]]
        rev_primer_pos <- coords[["rev_primer"]]

        
        # forward primer
        forw_primer <- toupper(str_trim(input$forw_primer))
        
        # reverse primer
        rev_primer <- toupper(str_trim(input$rev_primer))
        
        # coordinates of the first and last target codon nucleotide change
        firstCodonDifference <- min(codon_muts)
        lastCodonDifference <- max(codon_muts)
        
        # wild-type and mutant sequences making sure that the coordinates refer to them
        # "sequence" variable contains the mutant sequence
        
        # wild-type sequence
        # define the wild-type site assay for checking non-cutters
        wt_seq <- coords[["sequence"]]
        
        # subset the full sequence to that amplified by primers
        wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
        
        ## 2. perform design of the AS-PCR primers
        
        # add full names of all primers to the output based on the gene, mutation, direction and detection target 
        
        # add a list of mutated nucleotides within a primer to highlight them later during the output stage

        ################################################################################################################
        
        # run design functions
        forward_primers <-  design_forward_primers(sequence, wt_site_assay, lastCodonDifference, rev_primer)
        reverse_primers <-  design_reverse_primers(sequence, wt_site_assay, firstCodonDifference, forw_primer)
        
        primer_tables <- primerTablesOutput(forward_primers, reverse_primers)
        
        
        # add the final tags
        outputHTML <- paste(outputHTML, primer_tables, "</div>", sep = "")
        
        
        HTML(outputHTML)
        
      }) # end of lapply - outputList
      
      # CASE 2: input$mutatePAM == "yes" & input$REsites == "no" 
    } else if ( input$mutatePAM == "yes" & input$REsites == "no" ) {
      
      # the output for PAM/sgRNA only mutations
      isolate({ PAMonlyList <- PAM_mutations() })
      
      # get back the position data to the original state before off-setting
      codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
      pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
      
      # iterate each oligo design
      lapply(1:length(PAMonlyList), function(j) {
        
        # get sequence
        sequence <- toupper(PAMonlyList[[j]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        PAM_muts <- c()
        sgRNA_muts <- c()

        # get new replacement codon sequence
        new_codon <- substr(sequence, codon_pos[1], codon_pos[3])
        
        # codon mutations
        codon_muts <- PAMonlyList[[j]][["codon_diffs_coords"]]
        
        # get PAM mutations if they are available
        if(PAMonlyList[[j]][["PAM_mutant_codon"]] != "none"){
          
          PAM_muts <- PAMonlyList[[j]][["PAM_mut_codon_diffs"]]
          
        }
        
        # get sgRNA mutations
        if(PAMonlyList[[j]][["sgRNA_mutations"]]){
          
          sgRNA_muts <- PAMonlyList[[j]][["sgRNA_mut_codon_diffs"]] 
          
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
        
        # define the start and end of oligos for the purposes of correct output
        if(input$oriented == "sense"){
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[1] - 3 - input$leftArmLength
          oligoEnd <- pam_pos[1] - 3 + input$rightArmLength - 1
          
        }else{
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[3] + 3 - input$leftArmLength + 1
          oligoEnd <- pam_pos[3] + 3 + input$rightArmLength 
          
        } 
        
        # organize all mutated positions into a single vector
        all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        
        # add the strategy HTML if the program is at the beginning of the lapply loop
        if(j == 1){
          # initialize the output HTML
          outputHTML <- HTML(paste(strategyHTML, "<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j), "<br/>", sep = ""))          
          
        } else{
          
          # initialize the output HTML
          outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        }
        
        
        ############################################ formatOligo function ##############################################
        
        outputHTML <- paste(outputHTML,  formatOligo(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated), sep="")
        
        ###########################################################################################
        
        # get all information on the relevant restriction enzyme site
        if("enzymes" %in% names(PAMonlyList[[j]])){
          
          # iterate over all known enzymes for that sequence
          for(enzymeData in PAMonlyList[[j]][["enzymes"]]){
            
            # store the restriction site information
            enzyme <- enzymeData[["RE_enzyme"]]
            site <- enzymeData[["RE_site"]]
            site_coords <- enzymeData[["RE_site_coords"]]
            site_oriented <- enzymeData[["RE_site_oriented"]]
            
            
            # change the coordinates of the restriction site according to the oligo orientation
            if(input$orientedOligo == "anti"){
              site_coords <- nchar(sequence) - rev(site_coords) + 1
            }
        
            # offset all the positions
              
            # first, get the primers
            forw_primer_pos <- coords[["forw_primer"]]
            rev_primer_pos <- coords[["rev_primer"]]
            
            offset <- forw_primer_pos[1] -1
            
            # update coordinates and sequence to avoid incorrect calculation
            site_coords_offset <- enzymeData[["RE_site_coords"]] - offset
              
                
            # subset the full sequence to that amplified by primers
            SA_digestion <- substr(sequence, forw_primer_pos[1], rev_primer_pos[2])
              
            # add a second line for the restriction site
            outputHTML <- paste(outputHTML, "<span style='opacity: 0;'>", substr(sequence, oligoStart, site_coords[1]-1), "</span>", 
                                "<strong>",substr(sequence, site_coords[1], site_coords[2]), "</strong>", " - ", "<span style='color:blue'>",
                                enzyme, "</span>", " (", "<span style=' font-weight: bold'>",  site,  "</span>", ") ",
                                calculateFragments(RE_SITES, enzyme, site_coords_offset, site_oriented, SA_digestion),"<br/>", sep = "")
            
          }
          
          
        } # end of if statement for enzymes
        
        ############################################ AS-PCR primers output #############################################
        
        # steps for outputting tables of primers and their characteristics
        
        ## 1. collect variables for AS-PCR primer designs
        forw_primer_pos <- coords[["forw_primer"]]
        rev_primer_pos <- coords[["rev_primer"]]
        
        
        # forward primer
        forw_primer <- toupper(str_trim(input$forw_primer))
        
        # reverse primer
        rev_primer <- toupper(str_trim(input$rev_primer))
        
        # coordinates of the first and last target codon nucleotide change
        offset <- forw_primer_pos[1] - 1
        
        firstCodonDifference <- min(codon_muts) - offset
        lastCodonDifference <- max(codon_muts) - offset
        
        # wild-type and mutant sequences making sure that the coordinates refer to them
        # "sequence" variable contains the mutant sequence
        
        # wild-type sequence
        # define the wild-type site assay for checking non-cutters
        wt_seq <- coords[["sequence"]]
        
        # subset the full sequence to that amplified by primers
        wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
        
        ## 2. perform design of the AS-PCR primers
        
        # add full names of all primers to the output based on the gene, mutation, direction and detection target 
        
        # add a list of mutated nucleotides within a primer to highlight them later during the output stage
        
        ################################################################################################################
        
        # run design functions
        forward_primers <-  design_forward_primers(sequence, wt_site_assay, lastCodonDifference, rev_primer)
        reverse_primers <-  design_reverse_primers(sequence, wt_site_assay, firstCodonDifference, forw_primer)
        
        
        primer_tables <- primerTablesOutput(forward_primers, reverse_primers)
        
        # add the final tags
        outputHTML <- paste(outputHTML, primer_tables, "</div>", sep = "")
        
        
        HTML(outputHTML)
        
      }) # end of lapply - PAMonlyList
      
      
    } else if ( input$mutatePAM == "no" & input$REsites == "yes" ) {
      # CASE 3: input$mutatePAM == "no" & input$REsites == "yes"
      
      # list of silent REsite mutations designs
      isolate({ noPAM_REsitesList <- noPAMmuts_REsiteSilentMuts() })
      
      # iterate over the list items
      lapply(1:length(noPAM_REsitesList), function(j) {
        
        # get sequence
        sequence <- toupper(noPAM_REsitesList[[j]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        REsite_muts <- c()
        
        # get new replacement codon sequence
        new_codon <- substr(sequence, codon_pos[1], codon_pos[3])
        
        # codon mutations
        codon_muts <- noPAM_REsitesList[[j]][["codon_diffs_coords"]]
        
        # get the codon mutation that introduced a restriction site
        if( "enzymes" %in% names(noPAM_REsitesList[[j]])){
          
          # test if "RE_site_codon_diffs" in the relevant list
          if("RE_site_codon_diffs" %in% names(noPAM_REsitesList[[j]][["enzymes"]][[1]])){
            REsite_muts <- noPAM_REsitesList[[j]][["enzymes"]][[1]][["RE_site_codon_diffs"]]
          }  
        } 
        
        
        # combine all mutations
        mutated <- c(codon_muts, REsite_muts)
        
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
          
        }
        
        
        # define the start and end of oligos for the purposes of correct output
        if(input$oriented == "sense"){
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[1] - 3 - input$leftArmLength
          oligoEnd <- pam_pos[1] - 3 + input$rightArmLength - 1
          
        }else{
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[3] + 3 - input$leftArmLength + 1
          oligoEnd <- pam_pos[3] + 3 + input$rightArmLength 
          
        }
        
        # organize all mutated positions into a single vector
        all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        # consider checking whether start and end of the oligo are less and more than any labeled positions
        # in the sequence and then updating them accordingly
        # Also check that neither of the coordinates is negative or beyond the sequence length,
        # change them accordingly
        
        
        # add the strategy HTML if the program is at the beginning of the lapply loop
        if(j == 1){
          # initialize the output HTML
          outputHTML <- HTML(paste(strategyHTML, "<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        } else{
          
          # initialize the output HTML
          outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        }
        
        
        ############################################ formatOligo function ##############################################
        
        outputHTML <- paste(outputHTML,  formatOligo(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated), sep="")
        
        ###########################################################################################
        
        # get all information on the relevant restriction enzyme sites
        if("enzymes" %in% names(noPAM_REsitesList[[j]])){
          
          # iterate over all known enzymes for that sequence
          for(enzymeData in noPAM_REsitesList[[j]][["enzymes"]]){
            
            # store the restriction site information
            enzyme <- enzymeData[["RE_enzyme"]]
            site <- enzymeData[["RE_site"]]
            site_coords <- enzymeData[["RE_site_coords"]]
            site_oriented <- enzymeData[["RE_site_oriented"]]
            
            # change the coordinates of the restriction site according to the oligo orientation
            if(input$orientedOligo == "anti"){
              site_coords <- nchar(sequence) - rev(site_coords) + 1
            }
            
            # add a second line for the restriction site
            outputHTML <- paste(outputHTML, "<span style='opacity: 0;'>", substr(sequence, oligoStart, site_coords[1]-1), "</span>", 
                                "<strong>",substr(sequence, site_coords[1], site_coords[2]), "</strong>", " - ", "<span style='color:blue'>",
                                enzyme, "</span>", " (", "<span style=' font-weight: bold'>",  site,  "</span>", ") ",
                                calculateFragments(RE_SITES, enzyme, enzymeData[["RE_site_coords"]], site_oriented, sequence), "<br/>", sep = "")
            
          }
          
        }
        
        ############################################ AS-PCR primers output #############################################
        
        # steps for outputting tables of primers and their characteristics
        
        ## 1. collect variables for AS-PCR primer designs
        forw_primer_pos <- coords[["forw_primer"]]
        rev_primer_pos <- coords[["rev_primer"]]
        
        
        # forward primer
        forw_primer <- toupper(str_trim(input$forw_primer))
        
        # reverse primer
        rev_primer <- toupper(str_trim(input$rev_primer))
        
        # coordinates of the first and last target codon nucleotide change
        offset <- forw_primer_pos[1] - 1
        
        firstCodonDifference <- min(codon_muts)
        lastCodonDifference <- max(codon_muts)
        
        # wild-type and mutant sequences making sure that the coordinates refer to them
        # "sequence" variable contains the mutant sequence
        
        # wild-type sequence
        # define the wild-type site assay for checking non-cutters
        wt_seq <- coords[["sequence"]]
        
        # subset the full sequence to that amplified by primers
        wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
        
        ## 2. perform design of the AS-PCR primers
        
        # add full names of all primers to the output based on the gene, mutation, direction and detection target 
        
        # add a list of mutated nucleotides within a primer to highlight them later during the output stage
        
        ################################################################################################################
        
        # run design functions
        forward_primers <-  design_forward_primers(sequence, wt_site_assay, lastCodonDifference, rev_primer)
        reverse_primers <-  design_reverse_primers(sequence, wt_site_assay, firstCodonDifference, forw_primer)
        
        
        primer_tables <- primerTablesOutput(forward_primers, reverse_primers)
        
        # add the final tags
        outputHTML <- paste(outputHTML, primer_tables, "</div>", sep = "")
        
        
        HTML(outputHTML)
        
      }) # end of lapply - noPAM_REsitesList
      
    } else { # CASE 4: input$mutatePAM == "no" & input$REsites == "no"
      
      # codon mutations list
      isolate({ CodonMutsList <- codonMutations() })
      
      # get back the position data to the original state before off-setting
      codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
      pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
      
      # iterate over each codon mutation
      lapply(1:length(CodonMutsList), function(j){
        
        # get sequence
        sequence <- toupper(CodonMutsList[[j]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        
        # codon mutations are the only ones in this case
        mutated <- CodonMutsList[[j]][["codon_diffs_coords"]]
        
        # get new replacement codon sequence
        new_codon <- substr(sequence, codon_pos[1], codon_pos[3])
        
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
          
        }
        
        # define the start and end of oligos for the purposes of correct output
        if(input$oriented == "sense"){
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[1] - 3 - input$leftArmLength
          oligoEnd <- pam_pos[1] - 3 + input$rightArmLength - 1
          
        }else{
          
          # define the start and end of oligo coordinates
          oligoStart <- pam_pos[3] + 3 - input$leftArmLength + 1
          oligoEnd <- pam_pos[3] + 3 + input$rightArmLength 
          
        } 
        
        # organize all mutated positions into a single vector
        all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        
        # add the strategy HTML if the program is at the beginning of the lapply loop
        if(j == 1){
          # initialize the output HTML
          outputHTML <- HTML(paste(strategyHTML, "<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        } else{
          
          # initialize the output HTML
          outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        }
        
        
        ############################################ formatOligo function ##############################################
        
        outputHTML <- paste(outputHTML,  formatOligo(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated), sep="")
        
        ###########################################################################################
        
        # add the final tags
        outputHTML <- paste(outputHTML, "<br/>", sep = "")
        
        # get all information on the relevant restriction enzyme sites
        if("enzymes" %in% names(CodonMutsList[[j]])){
          
          # iterate over all known enzymes for that sequence
          for(enzymeData in CodonMutsList[[j]][["enzymes"]]){
            
            # store the restriction site information
            enzyme <- enzymeData[["RE_enzyme"]]
            site <- enzymeData[["RE_site"]]
            site_coords <- enzymeData[["RE_site_coords"]]
            site_oriented <- enzymeData[["RE_site_oriented"]]
            
            # change the coordinates of the restriction site according to the oligo orientation
            if(input$orientedOligo == "anti"){
              site_coords <- nchar(sequence) - rev(site_coords) + 1
            }
            
            # offset all the positions
            
            # first, get the primers
            forw_primer_pos <- coords[["forw_primer"]]
            rev_primer_pos <- coords[["rev_primer"]]
            
            offset <- forw_primer_pos[1] -1
            
            # update coordinates and sequence to avoid incorrect calculation
            site_coords_offset <- enzymeData[["RE_site_coords"]] - offset
            
            
            # subset the full sequence to that amplified by primers
            SA_digestion <- substr(sequence, forw_primer_pos[1], rev_primer_pos[2])
            
            # add a second line for the restriction site
            outputHTML <- paste(outputHTML, "<span style='opacity: 0;'>", substr(sequence, oligoStart, site_coords[1]-1), "</span>", 
                                "<strong>",substr(sequence, site_coords[1], site_coords[2]), "</strong>", " - ", "<span style='color:blue'>",
                                enzyme, "</span>", " (", "<span style=' font-weight: bold'>",  site,  "</span>", ") ",
                                calculateFragments(RE_SITES, enzyme, site_coords_offset, site_oriented, SA_digestion),"<br/>", sep = "")
            
          }
          
        }
        
        ############################################ AS-PCR primers output #############################################
        
        # steps for outputting tables of primers and their characteristics
        
        ## 1. collect variables for AS-PCR primer designs
        forw_primer_pos <- coords[["forw_primer"]]
        rev_primer_pos <- coords[["rev_primer"]]
        
        
        # forward primer
        forw_primer <- toupper(str_trim(input$forw_primer))
        
        # reverse primer
        rev_primer <- toupper(str_trim(input$rev_primer))
        
        # coordinates of the first and last target codon nucleotide change
        offset <- forw_primer_pos[1] - 1
        
        firstCodonDifference <- min(mutated) - offset
        lastCodonDifference <- max(mutated) - offset
        
        # wild-type and mutant sequences making sure that the coordinates refer to them
        # "sequence" variable contains the mutant sequence
        
        # wild-type sequence
        # define the wild-type site assay for checking non-cutters
        wt_seq <- coords[["sequence"]]
        
        # subset the full sequence to that amplified by primers
        wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
        
        ## 2. perform design of the AS-PCR primers
        
        # add full names of all primers to the output based on the gene, mutation, direction and detection target 
        
        # add a list of mutated nucleotides within a primer to highlight them later during the output stage
        
        ################################################################################################################
        
        # run design functions
        forward_primers <-  design_forward_primers(sequence, wt_site_assay, lastCodonDifference, rev_primer)
        reverse_primers <-  design_reverse_primers(sequence, wt_site_assay, firstCodonDifference, forw_primer)
        
        
        primer_tables <- primerTablesOutput(forward_primers, reverse_primers)
        
        # add the final tags
        outputHTML <- paste(outputHTML, primer_tables, "</div>", sep = "")
        
        HTML(outputHTML)
        
      }) # end of lapply loop
      
    } # end of if-else statement series
    
  } else{ # Cas12a PAMs
    
    # PAM/sgRNA mutations are not done, but sites are introduced
    # input$mutatePAM is not relevant so we do not test it
    
    if ( input$REsites == "yes" ) {
      # CASE 3: input$mutatePAM == "no" & input$REsites == "yes"
      
      # list of silent REsite mutations designs
      isolate({ noPAM_REsitesList <- noPAMmuts_REsiteSilentMuts() })
      
      # iterate over the list items
      lapply(1:length(noPAM_REsitesList), function(j) {
        
        # get sequence
        sequence <- toupper(noPAM_REsitesList[[j]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        REsite_muts <- c()
        
        # codon mutations
        codon_muts <- noPAM_REsitesList[[j]][["codon_diffs_coords"]]
        
        # get new replacement codon sequence
        new_codon <- substr(sequence, codon_pos[1], codon_pos[3])
        
        # get the codon mutation that introduced a restriction site
        if( "enzymes" %in% names(noPAM_REsitesList[[j]])){
          
          # test if "RE_site_codon_diffs" in the relevant list
          if("RE_site_codon_diffs" %in% names(noPAM_REsitesList[[j]][["enzymes"]][[1]])){
            REsite_muts <- noPAM_REsitesList[[j]][["enzymes"]][[1]][["RE_site_codon_diffs"]]
          }  
        } 
        
        
        # combine all mutations
        mutated <- c(codon_muts, REsite_muts)
        
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
          
        }
        
       # define the start and end of oligo coordinates
       # for Cas12a enzymes
      
        if(input$orientedOligo == "sense"){
          oligoStart <- pam_pos[1] - 18 - input$leftArmLength + 1
          oligoEnd <- pam_pos[1] - 18 + input$rightArmLength
          
        } else{
          
          oligoStart <- pam_pos[length(pam_pos)] + 18 - input$leftArmLength
          oligoEnd <- pam_pos[length(pam_pos)] + 18 + input$rightArmLength - 1
        } 
        
          
       # organize all mutated positions into a single vector
       all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        # consider checking whether start and end of the oligo are less and more than any labeled positions
        # in the sequence and then updating them accordingly
        # Also check that neither of the coordinates is negative or beyond the sequence length,
        # change them accordingly
        
        
        # add the strategy HTML if the program is at the beginning of the lapply loop
        if(j == 1){
          # initialize the output HTML
          outputHTML <- HTML(paste(strategyHTML, "<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        } else{
          
          # initialize the output HTML
          outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        }
        
        
       ############################################ formatOligo function ##############################################
       
       outputHTML <- paste(outputHTML,  formatOligo(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated), sep="")
       
       ###########################################################################################
       
        # get all information on the relevant restriction enzyme sites
        
        if("enzymes" %in% names(noPAM_REsitesList[[j]])){
          
          # iterate over all known enzymes for that sequence
          for(enzymeData in noPAM_REsitesList[[j]][["enzymes"]]){
            
            # store the restriction site information
            enzyme <- enzymeData[["RE_enzyme"]]
            site <- enzymeData[["RE_site"]]
            site_coords <- enzymeData[["RE_site_coords"]]
            site_oriented <- enzymeData[["RE_site_oriented"]]
            
            # change the coordinates of the restriction site according to the oligo orientation
            if(input$orientedOligo == "anti"){
              site_coords <- nchar(sequence) - rev(site_coords) + 1
            }
            
            
            # add a second line for the restriction site
            outputHTML <- paste(outputHTML, "<span style='opacity: 0;'>", substr(sequence, oligoStart, site_coords[1]-1), "</span>", 
                                "<strong>",substr(sequence, site_coords[1], site_coords[2]), "</strong>", " - ", "<span style='color:blue'>",
                                enzyme, "</span>", " (", "<span style=' font-weight: bold'>",  site,  "</span>", ") ",
                                calculateFragments(RE_SITES, enzyme, enzymeData[["RE_site_coords"]], site_oriented, sequence),
                                "<br/>", sep = "")
            
          }
          
        }
        
       ############################################ AS-PCR primers output #############################################
       
       # steps for outputting tables of primers and their characteristics
       
       ## 1. collect variables for AS-PCR primer designs
       forw_primer_pos <- coords[["forw_primer"]]
       rev_primer_pos <- coords[["rev_primer"]]
       
       
       # forward primer
       forw_primer <- toupper(str_trim(input$forw_primer))
       
       # reverse primer
       rev_primer <- toupper(str_trim(input$rev_primer))
       
       # coordinates of the first and last target codon nucleotide change
       offset <- forw_primer_pos[1] - 1
       
       firstCodonDifference <- min(codon_muts)
       lastCodonDifference <- max(codon_muts)
       
       # wild-type and mutant sequences making sure that the coordinates refer to them
       # "sequence" variable contains the mutant sequence
       
       # wild-type sequence
       # define the wild-type site assay for checking non-cutters
       wt_seq <- coords[["sequence"]]
       
       # subset the full sequence to that amplified by primers
       wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
       
       ## 2. perform design of the AS-PCR primers
       
       # add full names of all primers to the output based on the gene, mutation, direction and detection target 
       
       # add a list of mutated nucleotides within a primer to highlight them later during the output stage
       
       ################################################################################################################
       
       # run design functions
       forward_primers <-  design_forward_primers(sequence, wt_site_assay, lastCodonDifference, rev_primer)
       reverse_primers <-  design_reverse_primers(sequence, wt_site_assay, firstCodonDifference, forw_primer)
       
       
       primer_tables <- primerTablesOutput(forward_primers, reverse_primers)
       
       # add the final tags
       outputHTML <- paste(outputHTML, primer_tables, "</div>", sep = "")        
        
       HTML(outputHTML)
        
      }) # end of lapply - noPAM_REsitesList
      
    } else { # CASE 4: input$REsites == "no"
      
      # codon mutations list
      isolate({ CodonMutsList <- codonMutations() })
      
      # get back the position data to the original state before off-setting
      codon_pos <- coords[["codon"]][1]: coords[["codon"]][2]
      pam_pos <- coords[["PAM"]][1]: coords[["PAM"]][2]
      
      # iterate over each codon mutation
      lapply(1:length(CodonMutsList), function(j){
        
        # get sequence
        sequence <- toupper(CodonMutsList[[j]][["site_assay"]])
        
        # collect all mutated positions
        mutated <- c()
        
        # codon mutations are the only ones in this case
        mutated <- CodonMutsList[[j]][["codon_diffs_coords"]]
        
        # get new replacement codon sequence
        new_codon <- substr(sequence, codon_pos[1], codon_pos[3])
        
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
          
        }
        
        # define the start and end of oligos for the purposes of correct output
        # for Cas12a enzymes
        if(input$orientedOligo == "sense"){
          oligoStart <- pam_pos[1] - 18 - input$leftArmLength + 1
          oligoEnd <- pam_pos[1] - 18 + input$rightArmLength
          
        } else{
          
          oligoStart <- pam_pos[length(pam_pos)] + 18 - input$leftArmLength
          oligoEnd <- pam_pos[length(pam_pos)] + 18 + input$rightArmLength - 1
          
        } 
        
        
        # organize all mutated positions into a single vector
        all_special_pos <- sort(unique(c(codon_pos, pam_pos,mutated)))
        
        
        # add the strategy HTML if the program is at the beginning of the lapply loop
        if(j == 1){
          # initialize the output HTML
          outputHTML <- HTML(paste(strategyHTML, "<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        } else{
          
          # initialize the output HTML
          outputHTML <- HTML(paste("<div class='jumbotron', style='width: 100%; word-wrap:break-word; display:inline-block; font-family: Courier New; font-size: 14px;'>",
                                   oligoHeader(input$gene, input$Mutation, orig_codon, new_codon, j),  "<br/>", sep = ""))          
          
        }
        
        
        ############################################ formatOligo function ##############################################
        
        outputHTML <- paste(outputHTML,  formatOligo(sequence, oligoStart, oligoEnd, all_special_pos, codon_pos, pam_pos, mutated), sep="")
        
        ###########################################################################################
        
        # get all information on the relevant restriction enzyme sites
        if("enzymes" %in% names(CodonMutsList[[j]])){
          
          # iterate over all known enzymes for that sequence
          for(enzymeData in CodonMutsList[[j]][["enzymes"]]){
            
            # store the restriction site information
            enzyme <- enzymeData[["RE_enzyme"]]
            site <- enzymeData[["RE_site"]]
            site_coords <- enzymeData[["RE_site_coords"]]
            site_oriented <- enzymeData[["RE_site_oriented"]]
            
            # change the coordinates of the restriction site according to the oligo orientation
            if(input$orientedOligo == "anti"){
              site_coords <- nchar(sequence) - rev(site_coords) + 1
            }
            
            # offset all the positions
            
            # first, get the primers
            forw_primer_pos <- coords[["forw_primer"]]
            rev_primer_pos <- coords[["rev_primer"]]
            
            offset <- forw_primer_pos[1] -1
            
            # update coordinates and sequence to avoid incorrect calculation
            site_coords_offset <- enzymeData[["RE_site_coords"]] - offset
            
            
            # subset the full sequence to that amplified by primers
            SA_digestion <- substr(sequence, forw_primer_pos[1], rev_primer_pos[2])
            
            # add a second line for the restriction site
            outputHTML <- paste(outputHTML, "<span style='opacity: 0;'>", substr(sequence, oligoStart, site_coords[1]-1), "</span>", 
                                "<strong>",substr(sequence, site_coords[1], site_coords[2]), "</strong>", " - ", "<span style='color:blue'>",
                                enzyme, "</span>", " (", "<span style=' font-weight: bold'>",  site,  "</span>", ") ",
                                calculateFragments(RE_SITES, enzyme, site_coords_offset, site_oriented, SA_digestion),"<br/>", sep = "")
            
          }
          
        }
        
        ############################################ AS-PCR primers output #############################################
        
        # steps for outputting tables of primers and their characteristics
        
        ## 1. collect variables for AS-PCR primer designs
        forw_primer_pos <- coords[["forw_primer"]]
        rev_primer_pos <- coords[["rev_primer"]]
        
        
        # forward primer
        forw_primer <- toupper(str_trim(input$forw_primer))
        
        # reverse primer
        rev_primer <- toupper(str_trim(input$rev_primer))
        
        # coordinates of the first and last target codon nucleotide change
        offset <- forw_primer_pos[1] - 1
        
        firstCodonDifference <- min(mutated) - offset
        lastCodonDifference <- max(mutated) - offset
        
        # wild-type and mutant sequences making sure that the coordinates refer to them
        # "sequence" variable contains the mutant sequence
        
        # wild-type sequence
        # define the wild-type site assay for checking non-cutters
        wt_seq <- coords[["sequence"]]
        
        # subset the full sequence to that amplified by primers
        wt_site_assay <- substr(wt_seq, forw_primer_pos[1], rev_primer_pos[2])
        
        ## 2. perform design of the AS-PCR primers
        
        # add full names of all primers to the output based on the gene, mutation, direction and detection target 
        
        # add a list of mutated nucleotides within a primer to highlight them later during the output stage
        
        ################################################################################################################
        
        # run design functions
        forward_primers <-  design_forward_primers(sequence, wt_site_assay, lastCodonDifference, rev_primer)
        reverse_primers <-  design_reverse_primers(sequence, wt_site_assay, firstCodonDifference, forw_primer)
        
        
        primer_tables <- primerTablesOutput(forward_primers, reverse_primers)
        
        # add the final tags
        outputHTML <- paste(outputHTML, primer_tables, "</div>", sep = "")
        
        HTML(outputHTML)
        
      }) # end of lapply loop
      
    } # end of if-else statement series
    
  } # end of else part  
        
  }) # end of renderUI
    
 }) # end of observeEvent
  
  # handlers to ensure that the page with new oligo designs will show from top to bottom
  observeEvent(input$leftArmLength, {
    
    if(input$run >= 1){
      shinyjs::runjs("window.scrollTo(0,0);")  
    }
  })
  
  observeEvent(input$rightArmLength, {
    
    if(input$run >= 1){
      shinyjs::runjs("window.scrollTo(0,0);")  
    }
    
  })
  
  observeEvent(input$orientedOligo, {
    
    if(input$run >= 1){
      shinyjs::runjs("window.scrollTo(0,0);")  
    }
  
  })
  
  observeEvent(input$mutatePAM, {
  
    if(input$run >= 1){
      shinyjs::runjs("window.scrollTo(0,0);")  
    }
  
  })
  
  observeEvent(input$REsites, {
    
    if(input$run >= 1){
      shinyjs::runjs("window.scrollTo(0,0);")  
    }
  
  })
  
  
} # end of server

# Run the application 
shinyApp(ui = ui, server = server)
