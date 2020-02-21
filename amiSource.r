#r sourc file for ami functions






### treemaps modified from PR2
pr2.env <- new.env()

pr2.env$taxo_levels = c("Kingdom", "Phylum", "Class", "Order", 
    "Family", "Genus", "Species")
pr2.env$taxo_levels_number = length(pr2.env$taxo_levels)

pr2_treemap <- function(pr2, taxo_rank) {
    
    # Define the levels
    level1 = pr2.env$taxo_levels[taxo_rank]
    level2 = pr2.env$taxo_levels[taxo_rank + 3]
    # Group
    pr2_class <- pr2 %>% group_by_(level1, level2) %>% summarise(sequence_number = sum(abund))
    
    # Do a simple treemap
    treemap::treemap(pr2_class, index = c(level1, level2), vSize = "sequence_number", 
        title = paste(pel.names[i], sep=""), asp = 2, lowerbound.cex.labels = 0.2, fontsize.labels = 12, 
        palette = oceColorsJet(30), format.legend = list(scientific = FALSE, big.mark = " "))
}


ami_treemap<-function(pr2, taxo_rank, sample_type) {
    # utilised with kind acknowldegement of PR2 and the taxomap team 
    # Define the levels
    level1 = pr2.env$taxo_levels[taxo_rank]
    level2 = pr2.env$taxo_levels[taxo_rank + 2]
    # Group
    pr2_class <- pr2 %>% group_by_(level1, level2) %>% summarise(sequence_number = sum(abund))
    
    # Do a simple treemap
    treemap::treemap(pr2_class, index = c(level1, level2), vSize = "sequence_number", 
                      asp = 1, lowerbound.cex.labels = 0.2, 
                     palette = mycols200, vColor = level2, format.legend = list(scientific = FALSE, big.mark = " "),fontsize.labels=c(18,12),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                     fontcolor.labels=c("white","orange"),    # Color of labels
                     fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                     bg.labels=c("transparent"),              # Background color of labels
                     align.labels=list(
                         c("center", "center"), 
                         c("left", "top")
                     ),                                   # Where to place labels in the rectangle?
                     overlap.labels=1,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                     inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
                     aspRatio=3,
                     title=paste(sample_type),
                     lwds=(0.5)
                    
                     
    )
}

mycols108<-c("#3a43bb", "#74349a", "#234fa1", "#7d440c", "#8d353b", "#654387", "#a11d1b", "#9e1063", "#4c5814", "#84395b", 
"#764823", "#8e2286", "#7a3a7c", "#405b1b", "#014db6", "#605400", "#a90046", "#00622b", "#8b4000", "#405f00", 
"#ae0037", "#00633f", "#675827", "#b60072", "#bc1408", "#0f6e00", "#bb0062", "#5a6500", "#ae3d00", "#94505a", 
"#cb002d", "#016ca9", "#953ec5", "#a55800", "#df006a", "#00883f", "#867800", "#6a8200", "#d64908", "#ab6574", 
"#ca3cbf", "#99774b", "#9a60ea", "#009176", "#0092a9", "#ff2c6f", "#009b65", "#8c74ff", "#4887ff", "#d97100", 
"#ff4a6c", "#0297df", "#b09000", "#d468ef", "#ff50b0", "#ff6236", "#3d98ff", "#ff6170", "#85a06d", "#49b020", 
"#01a8d9", "#02afa0", "#ff715d", "#91a900", "#54afa3", "#d89300", "#4ab0b5", "#ff7987", "#cc9d00", "#ee9200", 
"#b8a771", "#00c164", "#d39c7b", "#00bac5", "#7fbe17", "#ff8f56", "#ff87d5", "#ff9571", "#ff82ff", "#ff978f", 
"#efa600", "#00ccab", "#ff9f48", "#e5ae00", "#c2a9ff", "#e69bff", "#ffa636", "#ff95f6", "#a2bcff", "#00cff1", 
"#ffa5c8", "#ffb094", "#a8d45f", "#dfc62c", "#e9c075", "#bbcd93", "#c7ce18", "#83db67", "#6cd6dd", "#ffb779", 
"#7ada97", "#5edcac", "#16e1a0", "#b8d074", "#f4b2e2", "#ebc23c", "#5edf7c", "#d4c960")

mycols310<-c("#374aa7", "#535602", "#1e5f11", "#704b0b", "#485827", "#a0193d", "#9b2628", "#922d4c", "#a3152d", "#8f3610", 
"#8a1d97", "#3f44b5", "#81411b", "#883a36", "#325c2e", "#584692", "#75406c", "#94216f", "#911f7f", "#813965", 
"#753499", "#684f02", "#8021a9", "#743a89", "#763d78", "#694e1e", "#5142aa", "#9f1554", "#793980", "#772ca8", 
"#4d4e7f", "#942565", "#872e7c", "#1750a3", "#275e24", "#005497", "#823e41", "#7b4432", "#0d6024", "#583bb6", 
"#863a4a", "#6d4377", "#653e9c", "#862a8c", "#9c280f", "#87375a", "#8b306a", "#75491f", "#92332a", "#883c2d", 
"#006135", "#ac0400", "#006302", "#a4007a", "#005a91", "#603bbf", "#006349", "#a70071", "#a82100", "#af003c", 
"#355988", "#af0051", "#a4008c", "#9a3900", "#9e069e", "#884658", "#b9001e", "#015ca8", "#b5005f", "#b5006f", 
"#665f00", "#793bc1", "#914d00", "#725d00", "#007019", "#00659c", "#0157d8", "#c30032", "#b70094", "#bc2e00", 
"#b03e00", "#726331", "#be008d", "#2a7400", "#a94a00", "#ce0041", "#706a00", "#a135bc", "#cb006e", "#437500", 
"#c7008a", "#00795b", "#007c14", "#9741c8", "#906200", "#405ee5", "#d12121", "#d5006a", "#875f8b", "#da005f", 
"#6a7500", "#99614a", "#0167ef", "#007e69", "#006dd9", "#be2cb2", "#006ae9", "#557c00", "#e00b3f", "#00833d", 
"#857000", "#6b5ee7", "#388400", "#e00080", "#b35e00", "#0278cb", "#d14502", "#6b7f00", "#9e6e00", "#00895e", 
"#ae6900", "#538700", "#ee115b", "#0082bc", "#e31e9b", "#027de4", "#0282cb", "#995de6", "#6c7baf", "#ed3040", 
"#898000", "#008e82", "#ef1a8b", "#00932d", "#018c9f", "#bf4ed4", "#5c8a60", "#a0784f", "#b16c81", "#01944b", 
"#cb6600", "#d841c1", "#f13d3c", "#4b79ff", "#b66f76", "#e43ab6", "#cc4fd3", "#fb2685", "#3e9a00", "#ff2b75", 
"#3484ff", "#019d4b", "#8874ff", "#0091d4", "#009f35", "#f1542c", "#0295c3", "#019d8f", "#ff4360", "#b88300", 
"#ff3aa1", "#be67f0", "#d87400", "#388eff", "#ff478f", "#ae8c00", "#00a1a6", "#ff535d", "#f848bd", "#bc7fa3", 
"#c58077", "#ff49a5", "#ff5c50", "#db63e8", "#999b00", "#01ab71", "#84a100", "#ff5f7c", "#7a8fff", "#129cff", 
"#fa6d27", "#ff6188", "#f5731b", "#9f87ff", "#00b158", "#00b247", "#ca77ff", "#ff62ae", "#ff6f54", "#be7fff", 
"#ff6e66", "#c89400", "#56b11a", "#ff6f7b", "#e68700", "#ff60d3", "#02b66d", "#ff6bbc", "#839bff", "#99a770", 
"#00b698", "#ff72a0", "#74b401", "#25bc41", "#00b0e6", "#bd8cff", "#ff75b4", "#f28b01", "#ff7d8d", "#ff8355", 
"#02b5df", "#ff8270", "#dd929e", "#00b9be", "#e57fff", "#d28aff", "#ff883c", "#00c07f", "#00b4f6", "#fb75f8", 
"#ff7ad9", "#00bcce", "#ff8e27", "#ff79e8", "#6dafff", "#ba9dff", "#b8a4da", "#5bc438", "#ff9061", "#00c59b", 
"#eaa000", "#ff9089", "#e39bb3", "#ff89c3", "#90bf0e", "#9dabff", "#00cd5e", "#87b3ff", "#ff93a3", "#cfb300", 
"#a3b3ea", "#ff9e40", "#ff9a8b", "#89c725", "#df9bff", "#ffa029", "#6abdff", "#02c8e6", "#ff9e78", "#d9b083", 
"#ff91f4", "#ffa35d", "#01d479", "#ff9cbb", "#adb4ff", "#bfafff", "#ffa83e", "#e5b500", "#ffa49b", "#ff9dd8", 
"#ffa4a6", "#f3a8ba", "#56d654", "#ff9bee", "#ffa1ce", "#fab00f", "#b9ca13", "#ffaf80", "#dfafff", "#66ceff", 
"#ffb45f", "#f8a8ff", "#ffb29a", "#c6cb78", "#7ad8a9", "#ebb2f6", "#5ed5ef", "#e7c32c", "#d6c76c", "#54de94", 
"#60d9c8", "#ffb93d", "#76dd5a", "#2ce18f", "#99d768", "#becd88", "#d2bbfe", "#2fdae5", "#cec88f", "#f6bb75", 
"#ffb3ab", "#d8c84f", "#6ddba1", "#36e277", "#a7d559", "#94d952", "#ebc15a", "#f6ba8a", "#61df73", "#ffade1", 
"#b5d32b", "#aad37c", "#b7d165", "#77dc7e", "#9dd678", "#fbadef", "#dbc583", "#e8c088", "#aad290", 
"#cfcc2f")
