In this tutorial we show how to set treatment contrasts in `<DESeq2>`
using the design or model matrix. This is a general and flexible way to
define contrasts, and is often useful for more complex contrasts or when
the design of the experiment is imbalanced (e.g. different number of
replicates in each group). Although we focus on `<DESeq2>`, the approach
can also be used with the other popular package `<edgeR>`.

Each section below covers a particular experimental design, from simpler
to more complex ones. The first chunk of code in each section is to
simulate data, which has no particular meaning and is only done in order
to have a DESeqDataSet object with the right kind of variables for each
example. In practice, users can ignore this step as they should have
created a DESeqDataSet object from their own data following the
[instructions in the
vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#the-deseqdataset).

There is a set of [accompanying
slides](https://docs.google.com/presentation/d/1B9zW1_F-kBqQEu4xqxIJrudYP5DecytYMRR6bY4H6aM/edit?usp=sharing)
that illustrate each of the sections below.

One factor, two levels (slide 5)
================================

    # simulate data
    dds <- makeExampleDESeqDataSet(n = 1000, m = 6, betaSD = 2)
    dds$condition <- factor(rep(c("shade", "sun"), each = 3))

First we can look at our sample information:

    colData(dds)

    ## DataFrame with 6 rows and 1 column
    ##         condition
    ##          <factor>
    ## sample1     shade
    ## sample2     shade
    ## sample3     shade
    ## sample4     sun  
    ## sample5     sun  
    ## sample6     sun

Our factor of interest is `condition` and so we define our design and
run the DESeq model fitting routine:

    design(dds) <- ~ 1 + condition # or just `~ condition`
    dds <- DESeq(dds) # equivalent to edgeR::glmFit()

Then check what coefficients DESeq estimated:

    resultsNames(dds)

    ## [1] "Intercept"              "condition_sun_vs_shade"

We can see that we have a coefficient for our *intercept* and
coefficient for the effect of “sun” (i.e. differences between sun versus
shade).

Using the more standard syntax, we can obtain the results for the effect
of sun as such:

    res1 <- results(dds, contrast = list("condition_sun_vs_shade"))
    res1

    ## log2 fold change (MLE): condition_sun_vs_shade effect 
    ## Wald test p-value: condition_sun_vs_shade effect 
    ## DataFrame with 1000 rows and 6 columns
    ##           baseMean log2FoldChange     lfcSE         stat      pvalue
    ##          <numeric>      <numeric> <numeric>    <numeric>   <numeric>
    ## gene1     40.90941    1.267525859  0.574144  2.207679752   0.0272666
    ## gene2     12.21876   -0.269917301  1.111127 -0.242922069   0.8080658
    ## gene3      1.91439   -3.538133611  2.564464 -1.379677442   0.1676860
    ## gene4     10.24472    0.954811627  1.166408  0.818591708   0.4130194
    ## gene5     13.16824    0.000656519  0.868780  0.000755679   0.9993971
    ## ...            ...            ...       ...          ...         ...
    ## gene996   40.43827     -1.0291276  0.554587    -1.855664 0.063501471
    ## gene997   52.88360      0.0622133  0.561981     0.110704 0.911851377
    ## gene998   73.06582      1.3271896  0.576695     2.301373 0.021370581
    ## gene999    8.87701     -5.8385374  1.549471    -3.768084 0.000164506
    ## gene1000  37.06533      1.2669314  0.602010     2.104501 0.035334764
    ##                 padj
    ##            <numeric>
    ## gene1      0.0712378
    ## gene2      0.8779871
    ## gene3      0.2943125
    ## gene4      0.5692485
    ## gene5      0.9996728
    ## ...              ...
    ## gene996  0.138827354
    ## gene997  0.948279388
    ## gene998  0.059599481
    ## gene999  0.000914882
    ## gene1000 0.087737235

The above is a simple way to obtain the results of interest. But it is
worth understanding how DESeq is getting to these results by looking at
the model’s matrix. DESeq defines the model matrix using base R
functionality:

    model.matrix(design(dds), colData(dds))

    ##         (Intercept) conditionsun
    ## sample1           1            0
    ## sample2           1            0
    ## sample3           1            0
    ## sample4           1            1
    ## sample5           1            1
    ## sample6           1            1
    ## attr(,"assign")
    ## [1] 0 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$condition
    ## [1] "contr.treatment"

We can see that R coded “condition” as a dummy variable, with an
intercept (common to all samples) and a “conditionsun” variable, which
adds the effect of sun to samples 4-6.

We can actually set our contrasts in `DESeq2::results()` using a numeric
vector. The way it works is to define a vector of “weights” for the
coefficient(s) we want to test for. In this case, we have `(Intercept)`
and `conditionsun` as our coefficients (see model matrix above), and we
want to test for the effect of sun, so our contrast vector would be
`c(0, 1)`. In other words, we don’t care about the value of
`(Intercept)` (so it has a weight of 0), and we’re only interested in
the effect of sun (so we give it a weight of 1).

In this case the design is very simple, so we could define our contrast
vector “manually”. But for complex designs this can get more difficult
to do, so it’s worth mentioning the general way in which we can define
this. For any contrast of interest, we can follow three steps:

-   Get the model matrix
-   Subset the matrix for each group of interest and calculate its
    column means - this results in a vector of coefficients for each
    group
-   Subtract the group vectors from each other according to the
    comparison we’re interested in

Let’s see this example in action:

    # get the model matrix
    mod_mat <- model.matrix(design(dds), colData(dds))
    mod_mat

    ##         (Intercept) conditionsun
    ## sample1           1            0
    ## sample2           1            0
    ## sample3           1            0
    ## sample4           1            1
    ## sample5           1            1
    ## sample6           1            1
    ## attr(,"assign")
    ## [1] 0 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$condition
    ## [1] "contr.treatment"

    # calculate the vector of coefficient weights in the sun
    sun <- colMeans(mod_mat[dds$condition == "sun", ])
    sun

    ##  (Intercept) conditionsun 
    ##            1            1

    # calculate the vector of coefficient weights in the shade
    shade <- colMeans(mod_mat[dds$condition == "shade", ])
    shade

    ##  (Intercept) conditionsun 
    ##            1            0

    # The contrast we are interested in is the difference between sun and shade
    sun - shade

    ##  (Intercept) conditionsun 
    ##            0            1

That last step is where we define our contrast vector, and we can give
this directly to the `results` function:

    # get the results for this contrast
    res2 <- results(dds, contrast = sun - shade)

This gives us exactly the same results as before, which we can check for
example by plotting the log-fold-changes between the first and second
approach:

    plot(res1$log2FoldChange, res2$log2FoldChange)

Extra: recoding the design (slide 12)
-------------------------------------

Often, we can use different model matrices that essentially correspond
to the same design. For example, we could recode our design above by
removing the intercept:

    design(dds) <- ~ 0 + condition
    dds <- DESeq(dds)
    resultsNames(dds)

    ## [1] "conditionshade" "conditionsun"

In this case we get a coefficient corresponding to the average
expression in shade and the average expression in the sun (rather than
the *difference* between sun and shade).

If we use the same contrast trick as before (using the model matrix), we
can see the result is the same:

    # get the model matrix
    mod_mat <- model.matrix(design(dds), colData(dds))
    mod_mat

    ##         conditionshade conditionsun
    ## sample1              1            0
    ## sample2              1            0
    ## sample3              1            0
    ## sample4              0            1
    ## sample5              0            1
    ## sample6              0            1
    ## attr(,"assign")
    ## [1] 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$condition
    ## [1] "contr.treatment"

    # calculate weights for coefficients in each condition
    sun <- colMeans(mod_mat[which(dds$condition == "sun"), ])
    shade <- colMeans(mod_mat[which(dds$condition == "shade"), ])

    # get the results for our contrast
    res3 <- results(dds, contrast = sun - shade)

Again, the results are essentially the same:

    plot(res1$log2FoldChange, res3$log2FoldChange)

In theory there’s no difference between these two ways of defining our
design. The design with an intercept is more common, but for the
purposes of understanding what’s going on, it’s sometimes easier to look
at models without intercept.

One factor, three levels (slide 6)
==================================

    # simulate data
    dds <- makeExampleDESeqDataSet(n = 1000, m = 9, betaSD = 2)
    dds$condition <- NULL
    dds$colour <- factor(rep(c("pink", "yellow", "white"), each = 3))
    dds$colour <- relevel(dds$colour, "white")

First we can look at our sample information:

    colData(dds)

    ## DataFrame with 9 rows and 1 column
    ##           colour
    ##         <factor>
    ## sample1   pink  
    ## sample2   pink  
    ## sample3   pink  
    ## sample4   yellow
    ## sample5   yellow
    ## sample6   yellow
    ## sample7   white 
    ## sample8   white 
    ## sample9   white

As in the previous example, we only have one factor of interest,
`condition`, and so we define our design and run the DESeq as before:

    design(dds) <- ~ 1 + colour
    dds <- DESeq(dds)

    # check the coefficients estimated by DEseq
    resultsNames(dds)

    ## [1] "Intercept"              "colour_pink_vs_white"   "colour_yellow_vs_white"

We see that now we have 3 coefficients:

-   “Intercept” corresponds to white colour (our reference level)
-   “colour\_pink\_vs\_white” corresponds to the difference between the
    reference level and pink
-   “colour\_yellow\_vs\_white” corresponds to the difference between
    the reference level and yellow

We could obtain the difference between white and any of the two colours
easily:

    res1_pink_white <- results(dds, contrast = list("colour_pink_vs_white"))
    res1_yellow_white <- results(dds, contrast = list("colour_yellow_vs_white"))

For comparing pink vs yellow, however, we need to compare two
coefficients with each other to check whether they are themselves
different (check the slide to see the illustration). This is how the
standard DESeq syntax would be:

    res1_pink_yellow <- results(dds, contrast = list("colour_pink_vs_white", 
                                                     "colour_yellow_vs_white"))

However, following our three steps detailed in the first section, we can
define our comparisons from the design matrix:

    # define the model matrix
    mod_mat <- model.matrix(design(dds), colData(dds))
    mod_mat

    ##         (Intercept) colourpink colouryellow
    ## sample1           1          1            0
    ## sample2           1          1            0
    ## sample3           1          1            0
    ## sample4           1          0            1
    ## sample5           1          0            1
    ## sample6           1          0            1
    ## sample7           1          0            0
    ## sample8           1          0            0
    ## sample9           1          0            0
    ## attr(,"assign")
    ## [1] 0 1 1
    ## attr(,"contrasts")
    ## attr(,"contrasts")$colour
    ## [1] "contr.treatment"

    # calculate coefficient vectors for each group
    pink <- colMeans(mod_mat[dds$colour == "pink", ])
    white <- colMeans(mod_mat[dds$colour == "white", ])
    yellow <- colMeans(mod_mat[dds$colour == "yellow", ])

And we can now define any contrasts we want:

    # obtain results for each pairwise contrast
    res2_pink_white <- results(dds, contrast = pink - white)
    res2_pink_yellow <- results(dds, contrast = pink - yellow)
    res2_yellow_white <- results(dds, contrast = yellow - white)

    # plot the results from the two approaches to check that they are identical
    plot(res1_pink_white$log2FoldChange, res2_pink_white$log2FoldChange)
    plot(res1_pink_yellow$log2FoldChange, res2_pink_yellow$log2FoldChange)
    plot(res1_yellow_white$log2FoldChange, res2_yellow_white$log2FoldChange)

With this approach, we could even define a more unusual contrast, for
example to find genes that differ between pigmented and non-pigmented
samples:

    # define vector of coefficients for pigmented samples
    pigmented <- colMeans(mod_mat[dds$colour %in% c("pink", "yellow"),])

    # Our contrast of interest is
    pigmented - white

    ##  (Intercept)   colourpink colouryellow 
    ##          0.0          0.5          0.5

Notice the contrast vector in this case assigns a “weight” of 0.5 to
each of `colourpink` and `colouryellow`. This is equivalent to saying
that we want to consider the average of pink and yellow expression. In
fact, we could have also defined our contrast vector like this:

    # average of pink and yellow minus white
    (pink + yellow)/2 - white

    ##  (Intercept)   colourpink colouryellow 
    ##          0.0          0.5          0.5

To obtain our results, we use the `results()` function as before:

    # get the results between pigmented and white
    res2_pigmented <- results(dds, contrast = pigmented - white)

Extra: why not define a new group in our design matrix?
-------------------------------------------------------

For this last example (pigmented vs white), we may have considered
creating a new variable in our column data:

    dds$pigmented <- factor(dds$colour %in% c("pink", "yellow"))
    colData(dds)

    ## DataFrame with 9 rows and 3 columns
    ##           colour sizeFactor pigmented
    ##         <factor>  <numeric>  <factor>
    ## sample1   pink     0.972928     TRUE 
    ## sample2   pink     0.985088     TRUE 
    ## sample3   pink     0.960749     TRUE 
    ## sample4   yellow   0.916582     TRUE 
    ## sample5   yellow   0.936918     TRUE 
    ## sample6   yellow   1.137368     TRUE 
    ## sample7   white    1.071972     FALSE
    ## sample8   white    1.141490     FALSE
    ## sample9   white    1.140135     FALSE

and then re-run DESeq with a new design:

    design(dds) <- ~ 1 + pigmented
    dds <- DESeq(dds)
    resultsNames(dds)

    ## [1] "Intercept"               "pigmented_TRUE_vs_FALSE"

    res1_pigmented <- results(dds, contrast = list("pigmented_TRUE_vs_FALSE"))

However, in this model the gene dispersion is estimated together for
pink and yellow samples as if they were replicates of each other, which
may result in inflated/deflated estimates. Instead, our approach above
estimates the error within each of those groups.

To check the difference one could compare the two approaches visually:

    # compare the log-fold-changes between the two approaches
    plot(res1_pigmented$log2FoldChange, res2_pigmented$log2FoldChange)
    abline(0, 1, col = "brown", lwd = 2)

    # compare the errors between the two approaches
    plot(res1_pigmented$lfcSE, res2_pigmented$lfcSE)
    abline(0, 1, col = "brown", lwd = 2)

Two factors with interaction (slide 7)
======================================

    # simulate data
    dds <- makeExampleDESeqDataSet(n = 1000, m = 12, betaSD = 2)
    dds$colour <- factor(rep(c("pink", "white"), each = 6))
    dds$colour <- relevel(dds$colour, "white")
    dds$condition <- factor(rep(c("sun", "shade"), 6))
    dds <- dds[, order(dds$colour, dds$condition)]
    colnames(dds) <- paste0("sample", 1:ncol(dds))

First let’s look at our sample information:

    colData(dds)

    ## DataFrame with 12 rows and 2 columns
    ##          condition   colour
    ##           <factor> <factor>
    ## sample1      shade    white
    ## sample2      shade    white
    ## sample3      shade    white
    ## sample4      sun      white
    ## sample5      sun      white
    ## ...            ...      ...
    ## sample8      shade     pink
    ## sample9      shade     pink
    ## sample10     sun       pink
    ## sample11     sun       pink
    ## sample12     sun       pink

This time we have two factors of interest, and we want to model both
with an interaction (i.e. we assume that white and pink samples may
respond differently to sun/shade). We define our design accordingly and
fit the model:

    design(dds) <- ~ 1 + colour + condition + colour:condition
    dds <- DESeq(dds)
    resultsNames(dds)

    ## [1] "Intercept"               "colour_pink_vs_white"   
    ## [3] "condition_sun_vs_shade"  "colourpink.conditionsun"

Because we have two factors and an interaction, the number of
comparisons we can do is larger. Using our three-step approach from the
model matrix, we do things exactly as we’ve been doing so far:

    # get the model matrix
    mod_mat <- model.matrix(design(dds), colData(dds))

    # Define coefficient vectors for each condition
    pink_shade <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "shade", ])
    pink_sun <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "sun", ])
    white_shade <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "shade", ])
    white_sun <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "sun", ])

We are now ready to define any contrast of interest from these vectors
(for completeness we show the equivalent syntax using the coefficient’s
names from DESeq).

Pink vs White (in the shade):

    res1 <- results(dds, contrast = pink_shade - white_shade)
    # or equivalently
    res2 <- results(dds, contrast = list("colour_pink_vs_white"))

Pink vs White (in the sun):

    res1 <- results(dds, contrast = pink_sun - white_sun)
    # or equivalently
    res2 <- results(dds, contrast = list(c("colour_pink_vs_white",
                                           "colourpink.conditionsun")))

Sun vs Shade (for whites):

    res1 <- results(dds, contrast = white_sun - white_shade)
    # or equivalently
    res2 <- results(dds, contrast = list(c("condition_sun_vs_shade")))

Sun vs Shade (for pinks):

    res1 <- results(dds, contrast = pink_sun - pink_shade)
    # or equivalently
    res2 <- results(dds, contrast = list(c("condition_sun_vs_shade", 
                                           "colourpink.conditionsun")))

Interaction between colour and condition (i.e. do pinks and whites
respond differently to the sun?):

    res1 <- results(dds, 
                    contrast = (pink_sun - pink_shade) - (white_sun - white_shade))
    # or equivalently
    res2 <- results(dds, contrast = list("colourpink.conditionsun"))

In conclusion, although we can define these contrasts using DESeq
coefficient names, it is somewhat more explicit (and perhaps intuitive?)
what it is we’re comparing using matrix-based contrasts.

Three factors, with nesting (slide 8)
=====================================

    # simulate data
    dds <- makeExampleDESeqDataSet(n = 1000, m = 24, betaSD = 2)
    dds$colour <- factor(rep(c("white", "pink"), each = 12))
    dds$colour <- relevel(dds$colour, "white")
    dds$species <- factor(rep(LETTERS[1:4], each = 6))
    dds$condition <- factor(rep(c("sun", "shade"), 12))
    dds <- dds[, order(dds$colour, dds$species, dds$condition)]
    colnames(dds) <- paste0("sample", 1:ncol(dds))

First let’s look at our sample information:

    colData(dds)

    ## DataFrame with 24 rows and 3 columns
    ##          condition   colour  species
    ##           <factor> <factor> <factor>
    ## sample1      shade    white        A
    ## sample2      shade    white        A
    ## sample3      shade    white        A
    ## sample4      sun      white        A
    ## sample5      sun      white        A
    ## ...            ...      ...      ...
    ## sample20     shade     pink        D
    ## sample21     shade     pink        D
    ## sample22     sun       pink        D
    ## sample23     sun       pink        D
    ## sample24     sun       pink        D

Now we have three factors, but species is *nested* within colour (i.e. a
species is either white or pink, it cannot be both). Therefore, colour
is a linear combination with species (or, another way to think about it
is that colour is redundant with species). Because of this, we will
define our design without including “colour”, although later we can
compare groups of species of the same colour with each other.

    design(dds) <- ~ 1 + species + condition + species:condition
    dds <- DESeq(dds)
    resultsNames(dds)

    ## [1] "Intercept"              "species_B_vs_A"         "species_C_vs_A"        
    ## [4] "species_D_vs_A"         "condition_sun_vs_shade" "speciesB.conditionsun" 
    ## [7] "speciesC.conditionsun"  "speciesD.conditionsun"

Now it’s harder to define contrasts between groups of species of the
same colour using DESeq’s coefficient names (although still possible).
But using the model matrix approach, we do it in exactly the same way we
have done so far!

Again, let’s define our groups from the model matrix:

    # get the model matrix
    mod_mat <- model.matrix(design(dds), colData(dds))

    # define coefficient vectors for each group
    pink_shade <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "shade", ])
    white_shade <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "shade", ])
    pink_sun <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "sun", ])
    white_sun <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "sun", ])

It’s worth looking at some of these vectors, to see that they are
composed of weighted coefficients from different species. For example,
for “pink” species, we have equal contribution from “speciesC” and
“speciesD”:

    pink_shade

    ##           (Intercept)              speciesB              speciesC 
    ##                   1.0                   0.0                   0.5 
    ##              speciesD          conditionsun speciesB:conditionsun 
    ##                   0.5                   0.0                   0.0 
    ## speciesC:conditionsun speciesD:conditionsun 
    ##                   0.0                   0.0

And so, when we define our contrasts, each species will be correctly
weighted:

    pink_sun - pink_shade

    ##           (Intercept)              speciesB              speciesC 
    ##                   0.0                   0.0                   0.0 
    ##              speciesD          conditionsun speciesB:conditionsun 
    ##                   0.0                   1.0                   0.0 
    ## speciesC:conditionsun speciesD:conditionsun 
    ##                   0.5                   0.5

We can set our contrasts in exactly the same way as we did in the
previous section (for completeness, we also give the contrasts using
DESeq’s named coefficients).

Pink vs White (in the shade):

    res1_pink_white_shade <- results(dds, contrast = pink_shade - white_shade)
    # or equivalently
    res2_pink_white_shade <- results(dds, 
                                     contrast = list(c("species_B_vs_A"),
                                                     c("species_C_vs_A",
                                                       "species_D_vs_A")))

Pink vs White (in the sun):

    res1_pink_white_sun <- results(dds, contrast = pink_sun - white_sun)
    # or equivalently
    res2_pink_white_sun <- results(dds, 
                               contrast = list(c("species_B_vs_A",
                                                 "speciesB.conditionsun"),
                                               c("species_C_vs_A", 
                                                 "species_D_vs_A",
                                                 "speciesC.conditionsun",
                                                 "speciesD.conditionsun")))

And so on, for other contrasts of interest…

Extra: imbalanced design
------------------------

Let’s take our previous example, but drop one of the samples from the
data, so that we only have 2 replicates for it.

    dds <- dds[, -13] # drop one of the species C samples
    dds <- DESeq(dds)
    resultsNames(dds)

    ## [1] "Intercept"              "species_B_vs_A"         "species_C_vs_A"        
    ## [4] "species_D_vs_A"         "condition_sun_vs_shade" "speciesB.conditionsun" 
    ## [7] "speciesC.conditionsun"  "speciesD.conditionsun"

Define our model matrix and coefficient vectors:

    mod_mat <- model.matrix(design(dds), colData(dds))
    mod_mat

    ##          (Intercept) speciesB speciesC speciesD conditionsun
    ## sample1            1        0        0        0            0
    ## sample2            1        0        0        0            0
    ## sample3            1        0        0        0            0
    ## sample4            1        0        0        0            1
    ## sample5            1        0        0        0            1
    ## sample6            1        0        0        0            1
    ## sample7            1        1        0        0            0
    ## sample8            1        1        0        0            0
    ## sample9            1        1        0        0            0
    ## sample10           1        1        0        0            1
    ## sample11           1        1        0        0            1
    ## sample12           1        1        0        0            1
    ## sample14           1        0        1        0            0
    ## sample15           1        0        1        0            0
    ## sample16           1        0        1        0            1
    ## sample17           1        0        1        0            1
    ## sample18           1        0        1        0            1
    ## sample19           1        0        0        1            0
    ## sample20           1        0        0        1            0
    ## sample21           1        0        0        1            0
    ## sample22           1        0        0        1            1
    ## sample23           1        0        0        1            1
    ## sample24           1        0        0        1            1
    ##          speciesB:conditionsun speciesC:conditionsun speciesD:conditionsun
    ## sample1                      0                     0                     0
    ## sample2                      0                     0                     0
    ## sample3                      0                     0                     0
    ## sample4                      0                     0                     0
    ## sample5                      0                     0                     0
    ## sample6                      0                     0                     0
    ## sample7                      0                     0                     0
    ## sample8                      0                     0                     0
    ## sample9                      0                     0                     0
    ## sample10                     1                     0                     0
    ## sample11                     1                     0                     0
    ## sample12                     1                     0                     0
    ## sample14                     0                     0                     0
    ## sample15                     0                     0                     0
    ## sample16                     0                     1                     0
    ## sample17                     0                     1                     0
    ## sample18                     0                     1                     0
    ## sample19                     0                     0                     0
    ## sample20                     0                     0                     0
    ## sample21                     0                     0                     0
    ## sample22                     0                     0                     1
    ## sample23                     0                     0                     1
    ## sample24                     0                     0                     1
    ## attr(,"assign")
    ## [1] 0 1 1 1 2 3 3 3
    ## attr(,"contrasts")
    ## attr(,"contrasts")$species
    ## [1] "contr.treatment"
    ## 
    ## attr(,"contrasts")$condition
    ## [1] "contr.treatment"

    # define coefficient vectors for each group
    pink_shade <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "shade", ])
    white_shade <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "shade", ])
    pink_sun <- colMeans(mod_mat[dds$colour == "pink" & dds$condition == "sun", ])
    white_sun <- colMeans(mod_mat[dds$colour == "white" & dds$condition == "sun", ])

Now let’s check what happens to the pink\_shade group:

    pink_shade

    ##           (Intercept)              speciesB              speciesC 
    ##                   1.0                   0.0                   0.4 
    ##              speciesD          conditionsun speciesB:conditionsun 
    ##                   0.6                   0.0                   0.0 
    ## speciesC:conditionsun speciesD:conditionsun 
    ##                   0.0                   0.0

Notice that whereas before “speciesC” and “speciesD” had each a weight
of 0.5, now they have different weights. That’s because for speciesC
there’s only 2 replicates. So, we have a total of 5 white individuals in
the shade (2 from species C and 3 from D). Therefore, when we calculate
the average coefficients for pinks, we need to do it as 0.4 x speciesC +
0.6 x speciesD.

The nice thing about this approach is that we do not need to worry about
any of this, the weights come from our `colMeans()` call automatically.
And now, any contrasts that we make will take these weights into
account:

    # pink vs white (in the shade)
    pink_shade - white_shade

    ##           (Intercept)              speciesB              speciesC 
    ##                   0.0                  -0.5                   0.4 
    ##              speciesD          conditionsun speciesB:conditionsun 
    ##                   0.6                   0.0                   0.0 
    ## speciesC:conditionsun speciesD:conditionsun 
    ##                   0.0                   0.0

    # interaction
    (pink_sun - pink_shade) - (white_sun - white_shade)

    ##           (Intercept)              speciesB              speciesC 
    ##                   0.0                   0.0                   0.1 
    ##              speciesD          conditionsun speciesB:conditionsun 
    ##                  -0.1                   0.0                  -0.5 
    ## speciesC:conditionsun speciesD:conditionsun 
    ##                   0.5                   0.5

Further reading
===============

-   Forum discussion about nested design:
    <a href="http://seqanswers.com/forums/showthread.php?t=47766" class="uri">http://seqanswers.com/forums/showthread.php?t=47766</a>
