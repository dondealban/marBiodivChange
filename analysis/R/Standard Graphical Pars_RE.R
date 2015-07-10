
set_graph_pars <- function(ptype = "panel4",pty="s")
{
    mgp <- c(2.5, 1, 0); mar <- c(4, 4, 2.5, 1); oma <- c(0, 0, 0, 0)
    switch(ptype,
           panel6 = par(mfrow=c(2,3), mgp=mgp, mar=mar, pty=pty,
                        oma=oma, bty="L", cex.lab =1.2),                  
           panel4 = par(mfrow=c(2,2), mgp=mgp, mar=mar, pty=pty,
                        oma=oma, bty="L", cex.lab =1.2),
           panel3 = par(mfrow=c(1,3), mgp=mgp, mar=mar, pty=pty,
                        oma=oma, bty="L", cex.axis=0.85),                        
           panel2 = par(mfrow=c(1,2), mgp=mgp, mar=mar, pty=pty,
                        oma=oma, bty="L", cex.axis=0.85),
           panel1 = par(mfrow=c(1,1), mgp=mgp, mar=mar, pty=pty,
                        oma=oma, bty="L", cex.axis=0.85))
}

add_panel_label <- function(ltype="a")
{
    text <- paste(LETTERS[letters==ltype])
    mtext(text=text, side=3, adj=0)
}

## change to 'TRUE' to look at the formatting...
if (FALSE) {

    plot_dummy <- function()
        plot(1:100, 1:100, type="n", xlab="x variable", ylab="y variable")
    ## 6 panels
    dev.new(width=8, height=7.5)
    set_graph_pars(ptype = "panel6")
    plot_dummy(); add_panel_label("a")
    plot_dummy(); add_panel_label("b")
    plot_dummy(); add_panel_label("c")
    plot_dummy(); add_panel_label("d")
    plot_dummy(); add_panel_label("e")
    plot_dummy(); add_panel_label("f")    
    ## 4 panels
    dev.new(width=8, height=7.5)
    set_graph_pars(ptype = "panel4")
    plot_dummy(); add_panel_label("a")
    plot_dummy(); add_panel_label("b")
    plot_dummy(); add_panel_label("c")
    plot_dummy(); add_panel_label("d")
    ## 3 panels
    dev.new(width=8, height=4)
    set_graph_pars(ptype = "panel3")
    plot_dummy(); add_panel_label("a")
    plot_dummy(); add_panel_label("b")
    plot_dummy(); add_panel_label("c")
    ## 2 panels
    dev.new(width=8, height=4)
    set_graph_pars(ptype = "panel2")
    plot_dummy(); add_panel_label("a")
    plot_dummy(); add_panel_label("b")
    ## 1 panel
    dev.new(width=8, height=4)
    set_graph_pars(ptype = "panel1")
    plot_dummy()
    ## clean up
    rm(plot_dummy)
}
