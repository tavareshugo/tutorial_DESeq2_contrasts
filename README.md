⚠️ **READ THIS BEFORE CONTINUING**

As I’ve come across more experimental designs, I have come to question the “universality” of the approach demonstrated in this tutorial. 
In particular, for partially crossed imbalanced designs (see [issue #9](https://github.com/tavareshugo/tutorial_DESeq2_contrasts/issues/9)), where this approach would in fact create a wrong contrast vector. 
When (if) I have time to revise this tutorial, I will come back to it. 
Until them, please see these instructions as personal notes from someone who is still learning and may get things wrong.

A general recommendation when working with more complex designs and custom contrasts is to always check the results output against the original data. 
For example, [plot the original data points](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts) for a given “significant” gene, to see if the expression in the individual samples makes sense with the `logFoldChange` reported from the model output. 
If there’s a big discrepancy, then it’s possible the contrast was not set correctly.

Despite the non-universality of this approach, I think the tutorial still gives a reasonable intuition of how contrast vectors work with DESeq2, at least for relatively simpler designs.

------------------------------------------------------------------------

This repository contains some code explaining how to set DESeq2 contrasts using the model matrix of our design. 

~~This approach is general and so works for a range of experimental designs, even more complex ones~~. (see warning above)

- [R code](DESeq2_contrasts.md)
- [slides](https://docs.google.com/presentation/d/1B9zW1_F-kBqQEu4xqxIJrudYP5DecytYMRR6bY4H6aM/edit?usp=sharing)

----

If something doesn't look correct, please let me know ([open an issue](https://github.com/tavareshugo/tutorial_DESeq2_contrasts/issues/new)). 
I would love to learn if I've gotten something wrong!

----

**Using these materials**

Feel free to use and adapt these materials to your own needs. 
If you do, please give attribution to: Hugo Tavares (University of Cambridge, Bioinformatics Training Facility)
