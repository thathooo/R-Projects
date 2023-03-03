# 1
pdf(file = "outputs/plots/1.1 netVisual.pdf", height = 8, width = 8)
netVisual_circle(cc_obj@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc_obj@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cc_obj, slot.name = "netP", measure = "weight")
dev.off()

# 2
pdf(file = "outputs/plots/1.2 netVisual_single.pdf", height = 12, width = 8)
par(mfrow = c(4,2), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# 3
pdf(file = "outputs/plots/1.3 netAnalysis_contribution.pdf", height = 8, width = 8)
netAnalysis_contribution(cc_obj, signaling = pwn)
dev.off()

# 4.1
pdf(file = "outputs/plots/1.4.1 netAnalysis_signalingRole_heatmap.pdf", height = 8, width = 8)
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "all")
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cc_obj, pattern = "incoming")
dev.off()

# # 4.2
# pdf(file = "outputs/plots/1.4.2 netAnalysis_signalingRole_heatmap.pdf", height = 8, width = 8)
# netAnalysis_signalingRole_scatter(cc_obj, signaling = pwn)
# netAnalysis_signalingRole_scatter(cc_obj, pattern = "outgoing")
# netAnalysis_signalingRole_scatter(cc_obj, pattern = "incoming")
# dev.off()

# 5
pdf(file = "outputs/plots/1.5 pathway_contribution.pdf", height = 6, width = 12)
netAnalysis_signalingRole_network(cc_obj, slot.name = "netP")
dev.off()

# 6
pdf(file = "outputs/plots/1.6.1 CommunicationPatterns_outgoing.pdf", height = 6, width = 12)
cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "outgoing", k = nPatterns)
# Visualize the communication pattern using river plot
netAnalysis_river(cc_obj, pattern = "outgoing")
# Visualize the communication pattern using dot plot
netAnalysis_dot(cc_obj, pattern = "outgoing")
dev.off()

# 6.2
pdf(file = "outputs/plots/1.6.2 CommunicationPatterns_incoming.pdf", height = 6, width = 12)
cc_obj <- identifyCommunicationPatterns(cc_obj, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cc_obj, pattern = "incoming")
# dot plot
netAnalysis_dot(cc_obj, pattern = "incoming")
dev.off()

# 7
