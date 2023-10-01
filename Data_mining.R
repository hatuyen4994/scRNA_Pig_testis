library(httr)
library(rvest)
library(XML)

##FROM ENTREZ to ENSEMBL ID
base.url  = "https://www.ncbi.nlm.nih.gov/gene/"
entrez_ID = 396660
url <- paste0(base.url, entrez_ID)
resp <- GET(url)
resp_body <- read_html(resp)
nodes = html_nodes(resp_body, xpath=".//a[@class='genome-browser-link']")
link = html_attr(nodes,  name = "href")
ensembl_ID = tail(strsplit(link,"/")[[1]],n=1)

##FROM ENSEMBL to ENTREZ

base.url  = "https://www.ncbi.nlm.nih.gov/gene/?term="
ensembl_ID = "ENSSSCG00010031912"
url <- paste0(base.url, ensembl_ID)
resp <- GET(url)
html <- read_html(resp)
text = html %>% html_elements('head') %>% html_element("title") %>% html_text()
gene_name = strsplit(text," ")[[1]][1]
Sys.sleep(3)



#ensembl_ID to orthologue genes of Sertoli
library(httr)
library(rvest)
library(XML)
to.convert = rowData(sce)$Symbol[grep("ENSSSCG",rowData(sce)$Symbol)]

base.url  = "https://www.ncbi.nlm.nih.gov/gene/?term="
gene.names = c()

for (ensembl_ID in to.convert) {
  url <- paste0(base.url, ensembl_ID)
  resp <- GET(url)
  html <- read_html(resp)
  text = html %>% html_elements('head') %>% html_element("title") %>% html_text()
  name = strsplit(text," ")[[1]][1]
  if (name=="No"){
    name=ensembl_ID    
  }
  gene.names = c(gene.names,name)
}
rowData(sce[to.convert])$Symbol = gene.names
rownames(sce) = rowData(sce)$Symbol
#df.out = data.frame(names = gene.names, row.names=to.convert)
#dim(df.out)
#write.table(df.out, file='Sert.ENS2Name.tsv', sep='\t', col.names = FALSE)

#ensembl_ID to orthologue genes
library(httr)
library(rvest)
library(XML)
to.convert = rowData(sce)$Symbol[grep("ENSSSCG",rowData(sce)$Symbol)]

base.url  = "https://www.ncbi.nlm.nih.gov/gene/?term="
gene.names = c()

for (ensembl_ID in to.convert) {
  url <- paste0(base.url, ensembl_ID)
  resp <- GET(url)
  html <- read_html(resp)
  text = html %>% html_elements('head') %>% html_element("title") %>% html_text()
  name = strsplit(text," ")[[1]][1]
  if (name=="No"){
    name=ensembl_ID    
  }
  gene.names = c(gene.names,name)
}
rowData(sce[to.convert])$Symbol = gene.names
rownames(sce) = rowData(sce)$Symbol
#df.out = data.frame(names = gene.names, row.names=to.convert)
#write.table(df.out, file='Ley.ENS2Name.tsv', sep='\t', col.names = FALSE)
