#Generating data
p=100;n=300;k=40;
mu = rep(0,p); Sigma = diag(p)
X = matrix(rnorm(n*p),n)
nonzero = 1:k
beta = 5*(1:p%in%nonzero)*sign(rnorm(p))/ sqrt(n)
y = X%*%beta + rnorm(n,1)

#Generate knockoff copy
Xk = create.gaussian(X,mu,Sigma)

#Generate importance statistic using knockoff package
W = stat.glmnet_coefdiff(X,Xk,y)

#Using filer_EM to obtain the final rejeciton set
U = 1:p #Use the location of the hypotheses as the side information
result = adaptiveKnockoff::filter_EM(W,U,df=1,alpha = c(0.1,0.2),reveal_prop = 0.3,mute = FALSE)
result = adaptiveKnockoff::filter_randomForest(W,z,alpha = c(0.1,0.2))

pp1 = adaptiveKnockoff::plot_ordering(res,nonzero = nonzero,start_index = 701)
pp1
pp2 = adaptiveKnockoff::plot_vanilla(W,alpha = seq(0.3,0.01,-0.01),nonzero = nonzero,start_index = 701)
prow <- plot_grid( pp1+ theme(legend.position="none"),
                   pp2 + theme(legend.position="none"),
                   #align = 'h',
                   #hjust = -1,
                   nrow = 2
)
plot_grid( prow, rel_heights = c(1,1))
