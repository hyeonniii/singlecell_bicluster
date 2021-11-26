# ADMM Algorithm Updates

## lagrangrian from
:
$
\begin{split}
        L_{\delta}(X,A,B,U,V) =
        & {1\over 2}\lVert Z-X \rVert^2_F +\gamma \lVert X_\Omega \rVert_F^2 +\sum_{i<j}P_{\lambda}(\lVert A_{ij}\rVert ) \\
        &+\sum_{k<l} P_{\lambda}(\lVert B_{kl} \rVert )+\sum_{i<j}U_{ij}^T(A_{ij}-( X_i-X_j ) ) + {\delta \over2} \sum_{i<j} \lVert{ A_{ij}-( X_i-X_j )}\rVert^2_2\\&+\sum_{k<l}V_{kl}^T(B_{kl}-( X^k-X^l ) ) + {\delta \over2}\sum_{k<l}\lVert{ B_{kl}-( X^k-X^l )}\rVert^2_2
\end{split}
$

## ADMM update

### update X 
$
\begin{split}
        X^{(k)}=
        \text{arg}\min_x {1 \over 2} 
        &\lVert{Z-X^{(k-1)}}\rVert^2_F+\gamma \lVert X^{(k-1)}_\Omega \rVert_F^2-\sum_{S}U_{ij}^T( {X_i - X_j} ) +{\delta\over2}\sum_{S} \lVert{ A_{ij}-( X_i-X_j )}\rVert^2_2\\
        &-\sum_{K}V_{ij}^T( {X^k- X^l})+{\delta \over2} \sum_{K} \lVert{ B_{kl}-( X^k-X^l )}\rVert^2_2
\end{split}$
\
$
\begin{split}
Vec(X)=\argmin_{vec(X)}{1 \over 2} 
        &\lVert{Vec(Z)-Vec(X^{(k-1)})}\rVert^2_F+\gamma \lVert E_{\Omega} Vec(X^{(k-1)}) \rVert_F^2-\sum_{S}U_{s}^TE_{ij}Vec(X^{(k-1)}) 
        \\& +{\delta\over2}\sum_{S} \lVert{ A_{ij}-E_{ij}Vec(X^{(k-1)})}\rVert^2_2-\sum_{N}V_{n}^TE_{kl}Vec(X^{(k-1)})\\&+{\delta \over2} \sum_{N} \lVert{ B_{ij}-E_{kl}Vec( X^{(k-1)})}\rVert^2_2
\end{split}
$ 

