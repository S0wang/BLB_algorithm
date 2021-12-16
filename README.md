# BLB_algorithm - solution for big data
The Bag of Little Bootstrap (BLB) Algorithms goes as follows:
\begin{enumerate}
    \item Take $m$ subsamples of size $b$ without replacement
    \begin{enumerate}
        \item $b$ can be a function of $n$; $b = n^\gamma; \gamma \in [0.5,1]$
        \item These subsamples could comprise a disjoint partition of the entire dataset. 
    \end{enumerate}
    \item For each subsample of size b, get 100 MC Estimations.
    \begin{enumerate}
        \item Create sample of weights = rowMeans(multinom(100, n, 1/b)).
        \item Calculate parameter of interest (theta) using weights.
    \end{enumerate}
    \item Take average of all calculated thetas. 

\end{enumerate}
