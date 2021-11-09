# Frequency Estimation

Some frequency estimation stuff in Matlab, mostly from the CRC for Robust and Adaptive Systems

This site presents some [Matlab (tm)](http://www.google.com/url?q=http%3A%2F%2Fwww.mathworks.com%2F&sa=D&sntz=1&usg=AFQjCNEN8sEnpi2M6VohHRH2Q4pK9OZEUA) code for estimation of the frequency of a single, constant tone in noise. While these methods may be extended to the multiharmonic and multi-tone cases, these programs do not include this extension.
Currently no explanations for the algorithms are available in HTML form. The only on-line reference is [here (PDF)](http://www.google.com/url?q=http%3A%2F%2Fespace.library.uq.edu.au%2Fview%2FUQ%3A10626&sa=D&sntz=1&usg=AFQjCNFGDc6N8gt32C3HAVY8n7mUlY8BCg). The reference list from this, with some linking of authors, is available here in HTML form.

Professor Barry Quinn has recently completed a book on this topic:

Quinn, B.G, Hannan, E.J. (2000) “[The Estimation and Tracking of Frequency](http://www.google.com/url?q=http%3A%2F%2Fwww.cambridge.org%2Fuk%2Fcatalogue%2Fcatalogue.asp%3Fisbn%3D0521804469&sa=D&sntz=1&usg=AFQjCNFIMBThtqEYttxzyJD3-31pxIoZ6Q)”

Barry has made matlab code for the algorithms he describes in the book available here.

| Algorithm Type                           | Matlab Code Links                        | Complexity   | Statistically Efficient? |
| ---------------------------------------- | ---------------------------------------- | ------------ | ------------------------ |
| Maximum Likelihood (ML) & Approximate ML | Maximum Likelihood                       | > O(T log T) | Yes                      |
|                                          | Periodogram Maximiser                    | > O(T log T) | Yes                      |
|                                          | Discrete-frequency Periodogram Maximizer | O(T log T)   | No                       |

| Fourier Coefficient Interpolation Techniques supplied by Eric Jacobsen, EF Data Corp (now Intel).
Some extra work done by Eric is available here (or mirrored here).

O(T log T) No
Signal Subspace Minimum Variance O(T3) No
Bartlett O(T3) No
Noise Subspace Pisarenko O(T) No
MUSIC O(T3) No
Phase Weighted Averagers
Lank-Reed-Pollon, Kay, Lovell-Williamson, Clarkson-Kootsookos-Quinn

O(T) Depends on technique. Most No.
Filtering Techniques Fernandes-Goodwin-de Souza O(T) Yes
Quinn-Fernandes O(T) Yes
Hannan-Huang ? ?
Nehorai-Porat ? ?
