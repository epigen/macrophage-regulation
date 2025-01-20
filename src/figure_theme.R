
### LIBRARIES
# select libraries
required_libs <- c(
    "ggplot2",
    "patchwork",
    "scales",
    "dplyr",
    "reshape2",
    "pheatmap",
    "ggrepel",
    "stringr",
    "gridGraphics",
    "ggplotify",
    "directlabels",
    "data.table",
    "svglite",
    "sna",
    "RColorBrewer",
     "Hmisc"
)

set.seed(42)

# load libraries
loaded_libs <- lapply(required_libs, function(x) suppressWarnings(suppressMessages(library(x, character.only = TRUE))))

options(stringsAsFactors=F)
                      
# default plotting theme
macro_theme <- function(){
    
    # settings
    font <- "Arial"
    size <- 8
    
    theme_bw(
        base_size=size,
        base_family = font
    ) %+replace% 
    
    theme(
      #grid elements
#       panel.grid.major = element_blank(),    #strip major gridlines
#       panel.grid.minor = element_blank(),    #strip minor gridlines
#       axis.ticks = element_blank(),          #strip axis ticks
      
#       strips axis lines ?
      
      #text elements
        text = element_text(              
                   family = font,           
                   size = size),
        
      plot.title = element_text(             #title
                   family = font,            #set font family
                   size = size,                #set font size
                   face = 'bold',            #bold typeface
                   hjust = 0.5,                #center align
                   vjust = 2),               #raise slightly
      
      plot.subtitle = element_text(          #subtitle
                   family = font,            #font family
                   size = size),               #font size
      
      plot.caption = element_text(           #caption
                   family = font,            #font family
                   size = size,                 #font size
                   hjust = 0.5),               #center align
      
      axis.title = element_text(             #axis titles
                   family = font,            #font family
                   size = size),               #font size
      
      axis.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = size),                #font size
        
        legend.text = element_text(              #axis text
                   family = font,            #axis famuly
                   size = size), 
        
        legend.title = element_text(              #axis text
                   family = font,            #axis famuly
                   size = size),
      
#       axis.text.x = element_text(            #margin for axis text
#                     margin=margin(5, b = 10))
    )
}

### FUNCTIONS

# function to break line after n=40 characters to shorten term names
addline_format <- function(x,...){
    return(gsub('(.{1,40})(\\s|\\.|$)', '\\1\n', x))
}
                      
# function to format treatmens
treatment_format <- function(x,...){
    x <- gsub('untreated', 'Untreated', x)
    
    x <- gsub('LO28', 'Listeria', x)
    
    x <- gsub('C_albicans', 'Candida', x) 
    x <- gsub('C albicans', 'Candida', x) 
    
    x <- gsub('IFN_gamma', 'IFN-g', x)
    x <- gsub('IFN gamma', "IFN-g", x)
    
    x <- gsub('IFN_beta', 'IFN-b', x)
    x <- gsub('IFN beta', 'IFN-b', x)
    
    x <- gsub('LCMV_Cl13', 'LCMV', x)
    x <- gsub('LCMV Cl13', 'LCMV', x)
    
    x <- gsub('_', ' ', x)
    x <- gsub('+', ' ', x, fixed=TRUE)
    x <- gsub('.', ' ', x, fixed=TRUE)
    
    x <- gsub('-2h', ' 2h', x)
    x <- gsub('-6h', ' 6h', x)
    x <- gsub('-24h', ' 24h', x)
    
    return(x)
}
                      
# function to format cluster numbers to romans
cluster_format <- function(x,...){
    x <- treatment_format(x)
    
    x <- gsub('cluster_10', 'X', x)
    x <- gsub('cluster 10', 'X', x)
    x <- gsub('10', 'X', x)
    
    x <- gsub('cluster_1', 'I', x)
    x <- gsub('cluster 1', 'I', x)
    x <- gsub('1', 'I', x)
    
    x <- gsub('cluster_2', 'II', x)
    x <- gsub('cluster 2', 'II', x)
    x <- gsub('2', 'II', x)
    
    x <- gsub('cluster_3', 'III', x)
    x <- gsub('cluster 3', 'III', x)
    x <- gsub('3', 'III', x)
    
    x <- gsub('cluster_4', 'IV', x)
    x <- gsub('cluster 4', 'IV', x)
     x <- gsub('4', 'IV', x)
    
    x <- gsub('cluster_5', 'V', x)
    x <- gsub('cluster 5', 'V', x)
    x <- gsub('5', 'V', x)
    
    x <- gsub('cluster_6', 'VI', x)
    x <- gsub('cluster 6', 'VI', x)
    x <- gsub('6', 'VI', x)
    
    x <- gsub('cluster_7', 'VII', x)
    x <- gsub('cluster 7', 'VII', x)
    x <- gsub('7', 'VII', x)
    
    x <- gsub('cluster_8', 'VIII', x)
    x <- gsub('cluster 8', 'VIII', x)
    x <- gsub('8', 'VIII', x)
    
    x <- gsub('cluster_9', 'IX', x)
    x <- gsub('cluster 9', 'IX', x)
    x <- gsub('9', 'IX', x)
    
    return(x)
}
                      
# function to format cluster numbers to A,B,C
cluster_format2 <- function(x,...){
    x <- treatment_format(x)
    
    x <- gsub('cluster_10', 'J', x)
    x <- gsub('cluster 10', 'J', x)
    x <- gsub('10', 'J', x)
    
    x <- gsub('cluster_1', 'A', x)
    x <- gsub('cluster 1', 'A', x)
    x <- gsub('1', 'A', x)
    
    x <- gsub('cluster_2', 'B', x)
    x <- gsub('cluster 2', 'B', x)
    x <- gsub('2', 'B', x)
    
    x <- gsub('cluster_3', 'C', x)
    x <- gsub('cluster 3', 'C', x)
    x <- gsub('3', 'C', x)
    
    x <- gsub('cluster_4', 'D', x)
    x <- gsub('cluster 4', 'D', x)
     x <- gsub('4', 'D', x)
    
    x <- gsub('cluster_5', 'E', x)
    x <- gsub('cluster 5', 'E', x)
    x <- gsub('5', 'E', x)
    
    x <- gsub('cluster_6', 'F', x)
    x <- gsub('cluster 6', 'F', x)
    x <- gsub('6', 'F', x)
    
    x <- gsub('cluster_7', 'G', x)
    x <- gsub('cluster 7', 'G', x)
    x <- gsub('7', 'G', x)
    
    x <- gsub('cluster_8', 'H', x)
    x <- gsub('cluster 8', 'H', x)
    x <- gsub('8', 'H', x)
    
    x <- gsub('cluster_9', 'I', x)
    x <- gsub('cluster 9', 'I', x)
    x <- gsub('9', 'I', x)
    
    return(x)
}
                      
remove_term_suffix <- function(db, terms){
    if(grepl('GO', db, fixed = TRUE)){
        # remove "(GO:.......)" from terms & abbreviate if necessary
        return(gsub("GO:.......", "", gsub("\\(GO:.......)", "", terms)))
    }
    if(grepl('WikiPathways', db, fixed = TRUE)){
        # remove WikiPathway IDs WP+numbers
        return(gsub("WP.*", "", terms))
    }
    
    return(terms)
}                  

                      
# extended ggsave
ggsave_new <- function(filename, results_path, plot, width=5, height=5, dpi=300){
    for (format in c('svg','png')){
        ggsave(
          paste0(filename,'.',format),
          plot = plot,
          device = format,
          path = file.path(results_path),
          scale = 1,
          dpi = dpi,
            width = width,
            height = height,
          limitsize = FALSE,
        )
    }
    
    ggsave(
          paste0(filename,'.pdf'),
          plot = plot,
          device = cairo_pdf,
          path = file.path(results_path),
          scale = 1,
          dpi = dpi,
            width = width,
            height = height,
          limitsize = FALSE,
        )
}
                      
                      
                      

                      
### DEFINITIONS (eg shape & color)
                      
## BULK
# set timepoint shapes
time_shapes <- c('0h' = 16,
             '2h' = 21,
             '4h' = 22,
             '6h' = 23,
             '8h' = 24,
             '24h' = 25
            )
                      
# set timepoint symbols as NUMBERS & SYMBOLS -> unused
time_shapes <- c('0h' = "\u2299",
             '2h' = "\u2461",
             '4h' = "\u2463",
             '6h' = "\u2465",
             '8h' = "\u2467",
             '24h' = "\u2297"
            )
# set timepoint sizes
time_sizes <- c('0h' = 0.5,
             '2h' = 1,
             '4h' = 2,
             '6h' = 3,
             '8h' = 4,
             '24h' = 5
            )
                      
# treatment colors
treatment_colors <- c('untreated'='#707070',
                      'C albicans'='#0072B2',
                      'IFN beta'='#D55E00',
                      'IFN gamma'='#E69F00',
                      'LCMV Cl13'='#CC79A7',
                      'LO28'='#009E73',
                      'LPS'='#56B4E9'
                      )                      
                      
# set timepoint colors
time_colors <- c('0h' = "#FFF7BC",
             '2h' = "#FFF7BC",
             '4h' = "#FEE391",
             '6h' = "#FEC44F",
             '8h' = "#FB9A29",
             '24h' = "#EC7014"
            )

# set cluster colors
cluster_colors <- c('1' = '#CC6677', 
                    '2' = '#332288', 
                    '3' = '#DDCC77', 
                    '4' = '#117733', 
                    '5' = '#88CCEE', 
                    '6' = '#882255', 
                    '7' = '#44AA99', 
                    '8' = '#999933', 
                    '9' = '#AA4499',
                    '10' = '#EE3377')
                      
cluster_alphas <- c('1' = 1, 
                    '2' = 0.9, 
                    '3' = 0.8, 
                    '4' = 0.7, 
                    '5' = 0.6, 
                    '6' = 0.5, 
                    '7' = 0.4, 
                    '8' = 0.3, 
                    '9' = 0.2,
                    '10' = 0.1)
                      
                      
# set colors for epigenetic correlation plots
correlation_colors <- c('RNA' = '#00FFFF',
                        'ATAC' = '#FF0000',
                        'correlated' = '#000000'
                       )
                      
                      
## SINGLE CELL
                      
# set condition colors
cond_colors <- c('untreated' = '#00BA38',
                 'LO28_2h' = '#619CFF',
                 'LO28_6h' = '#F8766D',
                 'LO28-6h' = '#F8766D',
                 'LO28-24h' = '#619CFF'
                )
                      
# small screen KO colors
KO_col <- c('mStat1'='#F8766D',
            'mIrf8'='#E68613',
            'mKdm1b'='#CD9600',
            'mNon-targeting'='#D3D3D3',
            'mHdac6'='#7CAE00',
            'mIrf9'='#0CB702',
            'mTyk2'='#00BE67',
            'mcJun'='#00C19A',
            'mStat2'='#00BFC4',
            'mEp300'='#00B8E7',
            'mKdm6b'='#00A9FF',
            'mCsf1r'='#8494FF',
            'mRela'='#C77CFF',
            'mSpi1'='#ED68ED',
            'mCreb1'='#FF61CC',
            'mJak1'='#FF68A1'
           )
                      
# set condition symbols
cond_shapes <- c('untreated' = 19,
             'LO28-6h' = 17,
             'LO28-24h' = 15
            )