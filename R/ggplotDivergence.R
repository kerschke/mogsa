#' @export
ggplotDivergence = function(grid) {
  positive.divergence = as.data.frame(grid$dec.space[which(grid$div > 0),])
  
  ggplotHeatmap(
    cbind.data.frame(grid$dec.space, grid$height),
    minimalistic.image = T
  ) + geom_tile(aes(x1, x2),
                data = positive.divergence
  )
}

# efficient = locallyNondominatedCPP(as.matrix(grid$height), grid$dims, T)
# quant = which(grid$height < quantile(grid$height, probs=c(0.05)))

# locally.efficient = as.data.frame(grid$dec.space[quant,])

# ggplot(asp = 1) +
#   geom_tile(aes(x1, x2, fill=I("red")),
#             data = negative.divergence
#   ) + geom_tile(aes(x1, x2, fill=I("black")),
#                 data = locally.efficient
#   ) + guides(fill = FALSE)
