#
# node.names <- data.frame(name = c(
#   paste0("N = ", format(nrow(dmps), big.mark=","), " DMPs"),
#   rownames(genic.df)))
#
# genic.df$source <- 0
# genic.df$target <- 1:nrow(genic.df)
#
# p.sankey <- sankeyNetwork(Links = genic.df,
#                           Nodes = node.names,
#                           Source = "source", Target = "target",
#                           Value = "value", NodeID = "name",
#                           fontSize = 12, nodeWidth = 20)
#
# f1 <- file.path(ODIR, "sankey.html")
# f2 <- file.path(ODIR, "sankey.png")
#
# saveNetwork(p.sankey, file = f1, selfcontained = F)
#
# # zoom increases the resolution it seems
# webshot::webshot(f1, f2, vwidth = 500, vheight = 250, zoom=5)
#
#
# # CpG island/shelf.shore --------------------------------------------------
#
# p.bar <- cpg.df %>%
#   filter(value >0.3) %>%
#   rownames_to_column("Location") %>%
#   dplyr::mutate(Percent = round(value)) %>%
#   ggplot(aes(fill = Location, x = 0, y = Percent)) +
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_manual(name = NULL, values = pal_locuszoom()(4)) +
#   coord_flip() +
#   theme_void() +
#   theme(plot.background = element_rect(fill = "transparent", colour = "transparent"),
#         panel.background = element_rect(fill = "transparent", colour = "transparent"),
#         plot.title = element_text(hjust = 0.5, size=20),
#         legend.text=element_text(size=16),
#         legend.position = "bottom",
#         plot.margin = margin(t=5, r=5,b=5,l=5)) +
#   guides(fill=guide_legend(ncol=2))
#
# cowplot::save_plot(plot = p.bar, filename=file.path(ODIR, "cpg_barchart.png"),
#                    base_width = 6, base_height = 3)
#
#
#
# # Pi chart ----------------------------------------------------------------
#
# N <- nrow(dmps)
#
# tally.df <- dmps %>%
#   dplyr::mutate(direction = ifelse(pi.diff < 0, "Hypo", "Hyper")) %>%
#   group_by(direction) %>%
#   summarize(prop = round(100 * n() / N)) %>%
#   arrange(-prop) %>%
#   dplyr::mutate(ypos = cumsum(prop) - 0.5*prop) %>%
#   dplyr::mutate(label = paste0(direction, "\n(", prop, "%)"))
#
#
# p.pi <- ggplot(tally.df, aes(x="", y=prop, fill=label)) +
#   geom_bar(stat="identity", width=1.5, color="white") +
#   coord_polar("y", start = 0, direction = -1) +
#   theme_void() +
#   theme(legend.position="none") +
#   geom_text(aes(y = ypos, label = label), color = "white", size = 4) +
#   scale_fill_manual(values = c(C.HYPER, C.HYPO))
#
# cowplot::save_plot(plot = p.pi, filename=file.path(ODIR, "pi.png"),
#                    base_width = 3, base_height = 3)
#
