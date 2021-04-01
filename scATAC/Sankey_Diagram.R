library(ArchR)
library(tibble)
library(dplyr)
library(tidyr)
library(networkD3)
library(htmlwidgets)
set.seed(1)

proj <- loadArchRProject(path = "HemeTutorial")
Cdata <- getCellColData(ArchRProj = proj, select = NULL, drop = FALSE)

label1 <- Cdata$superClusters
label2 <- as.character(Cdata$predictedGroup_un_E95)
label2[which(Cdata$superClusters=="colon" | Cdata$superClusters=="unidentified")] <- "colon&unidentified"
label2[which(label2=="pharynx")] <- "pharynx*"
names(label1) <- rownames(Cdata)
names(label2) <- rownames(Cdata)

#exclude colon and unidentified
celltypes1 <- c("pharynx","esophagus","lung","stomach","liver","pancreas","intestine")
celltypes2 <- c("pharynx*","esophagus-2","esophagus-1","respiratory-lung","respiratory-trachea","anaterior stomach","posterior stomach","liver-hepatocyte","hepatoblast","dorsal pancreas","duodenum")

links <- data.frame(source=0,target=0,value=0)
for(i in celltypes1){
    for(j in celltypes2){
        strength <- length(intersect(which(label1==i),which(label2==j)))
        if(strength==0){
            next
        }else{
            links <- rbind(links,c(i,j,strength))
        }
    }
}
links <- links[-1,]

nodes <- data.frame(name=c(celltypes1,celltypes2))

nodes$target=c(rep(TRUE,7),rep(FALSE,11))

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

pal=c("#E6C2DC","#F47D2B","#89C75F","#C06CAB","#90D5E4","#FEE500","#D51F26","#E6C2DC","#F47D2B","#D8A767","#89C75F","#208A42","#89288F","#C06CAB","#90D5E4","#8A9FD1","#FEE500","#D51F26")

my_color = 'd3.scaleOrdinal().domain([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]).range(["#E6C2DC","#F47D2B","#89C75F","#C06CAB","#90D5E4","#FEE500","#D51F26","#E6C2DC","#F47D2B","#D8A767","#89C75F","#208A42","#89288F","#C06CAB","#90D5E4","#8A9FD1","#FEE500","#D51F26"])'

# Make the Network
library(networkD3)
p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", iterations = 0,
              sinksRight=FALSE, colourScale = my_color, fontSize=14, nodeWidth = 15, height = 800, width = 800)

p$x$nodes$target <- nodes$target

p <- onRender(p,
  '
  function(el) {
    d3.select(el)
      .selectAll(".node text")
      .filter(d => d.target)
      .attr("x", -6)
      .attr("text-anchor", "end");
  } '
)

saveRDS(p,"Labels_sankey_diagram.rds")
