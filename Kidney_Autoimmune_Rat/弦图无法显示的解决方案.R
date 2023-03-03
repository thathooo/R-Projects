pdf(file ="figures/CADM_chord.pdf", width = 20, height =16)
for (i in 1:length(object.list)) {
  netVisual_aggregate(
    object.list[[i]], 
    signaling = "APP", 
    layout = "chord", 
    signaling.name = paste(
      pathways.show, 
      names(object.list)[i])
  )
}
dev.off()