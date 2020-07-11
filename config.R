# Notes to running rJava:
# Look for the error message when runing rJava where the library expects the Java file exits
# download the file from https://www.oracle.com/technetwork/java/javase/downloads/java-archive-javase11-5116896.html
# the missing file will then appear after installation here in this folder and then rJava can be run: /Library/Java/JavaVirtualMachines

# On linux
# sudo apt-get install -y default-jre
# sudo apt-get install -y default-jdk
# sudo R CMD javareconf
# install.packages("rJava")


# required packages
CRAN_packageList = c(
  "V8",
  "shiny",
  "readr",
  "tidyverse",
  "plotly",
  "testthat",
  "xlsx",
  "grid",
  "ggthemes",
  "openxlsx",
  "viridis",
  "scam",
  "shinydashboard",
  "shinyjs",
  "formattable",
  "plotly",
  "colourpicker",
  "scales",
  "extrafont",
  "shinyWidgets",
  "shinycssloaders",
  "RcppRoll",
  "shinyjqui",
  "forcats",
  "knitr",
  "bib2df",
  "DT",
  "vembedr",
  "MESS", #for auc
  "zip"
  )

biocPackageList <- c()

key_packages <- c("shiny","tidyverse")

# Custom palettes
pal_manual <- list(
  "pal_Cyril1" = c("#000080","#b03060","#006400","#00ff00","#ff00ff","#ff0000","#ffd700","#00bfff","#acb54c","#b57b4c","#4cb5af","#724cb5"),
  "pal_Cyril2" = c("#4363d8","#3cb44b","#e6194b","#911eb4","#9a6324","#f58231","#000075","#808000","#ffe119","#bcf60c","#46f0f0","#f032e6"),
  "pal_subway" = c("#EE473E","#A9BD38","#0172BE","#C71F65","#A26337","#03B498","#8D7EB9","#D0A95F","#04B0EF","#F7A527","#C7BFB2","#1AB268"),
  "colorblind" = c("#8F0104","#0368D3","#0A4C4B","#380686","#874B10","#8ACEC8","#BF66F5","#E06A11","#BF89A0","#73B2FC","#2DF628","#EBBAD7","#AEDCFA","#F9FF93"),
  "deuteranopia" = c("#C8433A","#4B7AB3","#2F5D53","#9B251C","#1F4C92","#213C39","#7BA73D","#CD7EBD","#64A996","#E1C59A","#84E9E0","#DDB7C2","#D4F7E1","#D7FEFA","#FAF3F5"),
  "protanopia" = c("#CF402B","#697BBC","#2F635F","#273517","#374992","#6B212F","#ED6C52","#64A6C0","#58A296","#EFAB98","#D9D8E7","#78D2C3","#F9E8DB","#F2F0FB","#B2FEF3"),
  "tritanopia" = c("#B02C54","#1D5399","#DA5E41","#3D8C45","#5E5782","#EB6545","#79ABD3","#977EDD","#EEAE99","#BFEE74","#C9C5E2","#F4FA42","#F1F5FE","#EEF4CE")
)


# Notes on color
Menubar_topleft = "#3C8B52"
Menubar_centerright = "#4AA361"
FullBoxheader = "#4AA361"
Slimboxheader = "#5E9676"
Page_background = "#EDF0F5"
Navbar_background = "#242D32"
valid_color_grey_boxborder = "#CACED7"

# Spinner
spinner_col = Menubar_centerright
spinner_bg = Page_background

# On Off switches
on_color = Menubar_centerright
off_color = "#8c4747"

# external links
beginners_tutorial            <- "https://youtu.be/nOAVuvtYCnk"
advanced_options              <- "https://youtu.be/nOAVuvtYCnk"
multiple_mixed_experiments    <- "https://youtu.be/nOAVuvtYCnk"
interpreting_results          <- "https://youtu.be/nOAVuvtYCnk"
R_interface                   <- "https://youtu.be/nOAVuvtYCnk"