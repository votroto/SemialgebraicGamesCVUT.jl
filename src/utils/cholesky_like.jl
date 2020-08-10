# Improving this algorithm will likely improve the performance of the
# Change-of-Basis model

""" Implements Cholesky-like decomposition X=VV' for PSD matrices """
function cholesky_like(X)
        u, s, v = svd(X)
        q, r = qr(sqrt(Diagonal(s)) * u')
        r
end