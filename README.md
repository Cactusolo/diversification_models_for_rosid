# diversification_models_for_rosid

**why the maximum-likelihood value is not negative?**
In Rabosky and Lovette (2008), the authors metioned the likelihood function (eq. 6) frequently resulted in positive log-likelihood values; this occurs because the (i − 1) lamda (t) P(ti, T) component of the likelihood is not normalized and is merely proportional to the actual probability density. 

And moreover, they could not use the likelihood-ratio test to compare models because the SPVAR and EXVAR models are not nested; rather, they compared model fits using the Akaike Information Criterion (AIC), and the lowest ΔAIC indicating the best-fit model.

**Akaike Weight**

Since it's difficult to unambiguously interpret the observed AIC differences, we transformed the “raw” AIC values to so-called Akaike weights (See Wagenmakers and Farrell, 2004). The Akaike weight is the more close to 1, the correspnded model is the more preferred model over all the models compared.

_Citation:_

Wagenmakers EJ, Farrell S 2004. AIC model selection using Akaike weights. Psychonomic Bulletin & Review 11: 192–196.

Rabosky DL, Lovette IJ. 2008. Explosive evolutionary radiations: decreasing speciation or increasing extinction through time? Evolution 62: 1866–1875.
