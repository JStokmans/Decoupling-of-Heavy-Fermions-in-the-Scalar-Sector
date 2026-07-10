* This code computes the expressions for each term in the Forest Formula for the Beetle Diagram and puts these over a common denominator so they can be added together later.
* It stores the resulting expressions in a file "Beetle_Diagram_Forests.sav", which is placed in the same directory as this code.
* Note that if the computation takes an unreasonably long time, one can choose to comment out several of the six permutations of the Beetle Diagram. 

Dimension 4;
Symbols a,n,k,m,M;
Symbols [lf^2-M^2],[(lf+p)^2-M^2],
        [(lf+l1)^2-M^2],[(lf+l2)^2-M^2],
        [(lf+p-l1)^2-M^2],[(lf+p-l2)^2-M^2],
        [(lf+p-l1+l2)^2-M^2],[(lf+l1-l2)^2-M^2],
        [(l2-l1)^2-m^2],[(l1-p)^2-m^2],[l2^2-m^2];      * Denominators of I_Gamma
Symbols [(lf-l2)^2-M^2];                                * Denominators of t_gamma_3*I_Gamma
Symbols [(Qf+l1)^2-M^2],[(Qf+l2)^2-M^2],
        [(Qf+Q-l1)^2-M^2],[(Qf+Q-l2)^2-M^2],
        [(Qf+Q-l1+l2)^2-M^2],[(Qf+l1-l2)^2-M^2],
        [(l1-Q)^2-m^2];                                 * Denominators of t_gamma_1*I_Gamma
Symbols [l1^2-m^2],
        [(lf-l1)^2-M^2],
        [(lf-l1+l2)^2-M^2];                             * Denominators of t_Gamma*\overline{R}_Gamma
Functions [d/dp],f (Symmetric);
CFunctions g (Symmetric);
Vectors x,y,z;
Vectors lf,l1,l2,p;         * Unfixed momenta
Vectors Q,Qf;               * gamma_1 renormalization scale
Indices i,j;

OFF Statistics;

#procedure FactorDenominators
Bracket
    [lf^2-M^2],
    [(lf+p)^2-M^2],
    [(lf+l1)^2-M^2],
    [(lf+l2)^2-M^2],
    [(lf+p-l1)^2-M^2],
    [(lf+p-l2)^2-M^2],
    [(lf+p-l1+l2)^2-M^2],
    [(lf+l1-l2)^2-M^2],
    [(l1-p)^2-m^2],
    [(l2-l1)^2-m^2],
    [l2^2-m^2],
    [(lf-l2)^2-M^2],
    [(Qf+l1)^2-M^2],
    [(Qf+l2)^2-M^2],
    [(Qf+Q-l1)^2-M^2],
    [(Qf+Q-l2)^2-M^2],
    [(Qf+Q-l1+l2)^2-M^2],
    [(Qf+l1-l2)^2-M^2],
    [(l1-Q)^2-m^2],
    [l1^2-m^2],
    [(lf-l1)^2-M^2],
    [(lf-l1+l2)^2-M^2]
    ;
#endprocedure

#procedure ReplaceMomenta
* Replacing mometa by noncommuting functions so they don't commute with the derivative operator.
    id x?.y? = f(x,y);
    id [(lf+p)^2-M^2]^n?        = f([(lf+p)^2-M^2]^n);
    id [(l1-p)^2-m^2]^n?        = f([(l1-p)^2-m^2]^n);
    id [(lf+p-l1)^2-M^2]^n?     = f([(lf+p-l1)^2-M^2]^n);
    id [(lf+p-l2)^2-M^2]^n?     = f([(lf+p-l2)^2-M^2]^n);
    id [(lf+p-l1+l2)^2-M^2]^n?  = f([(lf+p-l1+l2)^2-M^2]^n);
    id f(1) = 1;
    id x?(i?) = g(x(i));
#endprocedure

#procedure Derivative
Repeat;
* All occurrences of dot products.
    id g(x?(i?))*[d/dp](x?)*f(y?!{x?},z?!{x?})       = f(y,z)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(z?!{x?},g(y?))         = f(g(y),z)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(g(y?),g(z?))           = f(g(y),g(z))*g(x(i))*[d/dp](x);

    id g(x?(i?))*[d/dp](x?)*f(x?,x?)                 = 2*f(g(x),x)+f(x,x)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(y?!{x?},x?)            = g(g(x),y)+f(y,x)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(a?,x?)                 = g(g(x),a)+f(a,x)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(g(y?),x?)              = g(g(x),g(y))+f(g(y),x)*g(x(i))*[d/dp](x);

* All denominators.
    id g(x?(i?))*[d/dp](x?)*f([(lf+p)^2-M^2]^n?)        = (+1*n*2*(g(lf,g(x))+f(p,g(x)))*f([(lf+p)^2-M^2]^(n-1)) + f([(lf+p)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?(i?))*[d/dp](x?)*f([(l1-p)^2-m^2]^n?)        = (-1*n*2*(g(l1,g(x))-f(p,g(x)))*f([(l1-p)^2-m^2]^(n-1)) + f([(l1-p)^2-m^2]^n)*g(x(i))*[d/dp](x));
    id g(x?(i?))*[d/dp](x?)*f([(lf+p-l1)^2-M^2]^n?)     = (+1*n*2*(g(lf,g(x))+f(p,g(x))-g(l1,g(x)))*f([(lf+p-l1)^2-M^2]^(n-1)) + f([(lf+p-l1)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?(i?))*[d/dp](x?)*f([(lf+p-l2)^2-M^2]^n?)     = (+1*n*2*(g(lf,g(x))+f(p,g(x))-g(l2,g(x)))*f([(lf+p-l2)^2-M^2]^(n-1)) + f([(lf+p-l2)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?(i?))*[d/dp](x?)*f([(lf+p-l1+l2)^2-M^2]^n?)  = (+1*n*2*(g(lf,g(x))+f(p,g(x))-g(l1,g(x))+g(l2,g(x)))*f([(lf+p-l1+l2)^2-M^2]^(n-1)) + f([(lf+p-l1+l2)^2-M^2]^n)*g(x(i))*[d/dp](x));
EndRepeat;
id [d/dp](x?) = 0;
#endprocedure

#procedure RestoreMomenta
* Restoring the momenta, i.e. removing the functions f and g used to facilitate differentiation.
    id f(x?,y?) = x.y;
    id f(x?,g(y?)) = x.y;
    id f(g(x?),g(y?)) = x.y;
    id g(x?,y?) = x.y;
    id g(x?,g(y?)) = x.y;
    id g(g(x?),g(y?)) = x.y;

    id f(0,y?) = 0;
    id f(0,g(y?)) = 0;
    id f(g(0),g(y?)) = 0;
    id g(0,y?) = 0;
    id g(0,g(y?)) = 0;
    id g(g(0),g(y?)) = 0;

    id f(0,0) = 0;

    id f(a?) = a;
#endprocedure

* Defining the proper forests.
L [I_Gamma] =   (g_(1,lf)+M*gi_(1))*(

(g_(1,lf)+g_(1,l2)+M*gi_(1))*(g_(1,lf)+g_(1,l1)+M*gi_(1))/([(lf+l1)^2-M^2]*[(lf+l2)^2-M^2])
+(g_(1,lf)+g_(1,l1)-g_(1,l2)+M*gi_(1))*(g_(1,lf)+g_(1,l1)+M*gi_(1))/([(lf+l1-l2)^2-M^2]*[(lf+l1)^2-M^2])
+(g_(1,lf)+g_(1,l1)-g_(1,l2)+M*gi_(1))*(g_(1,lf)+g_(1,p)-g_(1,l2)+M*gi_(1))/([(lf+l1-l2)^2-M^2]*[(lf+p-l2)^2-M^2])
+(g_(1,lf)+g_(1,l2)+M*gi_(1))*(g_(1,lf)+g_(1,p)-g_(1,l1)+g_(1,l2)+M*gi_(1))/([(lf+l2)^2-M^2]*[(lf+p-l1+l2)^2-M^2])
+(g_(1,lf)+g_(1,p)-g_(1,l1)+M*gi_(1))*(g_(1,lf)+g_(1,p)-g_(1,l2)+M*gi_(1))/([(lf+p-l1)^2-M^2]*[(lf+p-l2)^2-M^2])
+(g_(1,lf)+g_(1,p)-g_(1,l1)+M*gi_(1))*(g_(1,lf)+g_(1,p)-g_(1,l1)+g_(1,l2)+M*gi_(1))/([(lf+p-l1)^2-M^2]*[(lf+p-l1+l2)^2-M^2])

                )*(g_(1,lf)+g_(1,p)+M*gi_(1)) / 

                ([lf^2-M^2]*[(lf+p)^2-M^2]*[(l2-l1)^2-m^2]*[(l1-p)^2-m^2]*[l2^2-m^2]);

L [t_gamma_f*I_Gamma] = (g_(1,lf)+M*gi_(1))*(

6*
(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])

                        )*(g_(1,lf)+M*gi_(1)) / 

                        ([lf^2-M^2]*[lf^2-M^2]*[(l2-l1)^2-m^2]*[(l1-p)^2-m^2]*[l2^2-m^2]);

L [t_gamma_3*I_Gamma] = (g_(1,lf)+M*gi_(1))*(

(g_(1,lf)+g_(1,l2)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[(lf+l2)^2-M^2])
+(g_(1,lf)-g_(1,l2)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([(lf-l2)^2-M^2]*[lf^2-M^2])
+(g_(1,lf)-g_(1,l2)+M*gi_(1))*(g_(1,lf)-g_(1,l2)+M*gi_(1))/([(lf-l2)^2-M^2]*[(lf-l2)^2-M^2])
+(g_(1,lf)+g_(1,l2)+M*gi_(1))*(g_(1,lf)+g_(1,l2)+M*gi_(1))/([(lf+l2)^2-M^2]*[(lf+l2)^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)-g_(1,l2)+M*gi_(1))/([lf^2-M^2]*[(lf-l2)^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)+g_(1,l2)+M*gi_(1))/([lf^2-M^2]*[(lf+l2)^2-M^2])

                        )*(g_(1,lf)+M*gi_(1)) / 

                        ([lf^2-M^2]*[lf^2-M^2]*[l2^2-m^2]*[l2^2-m^2]) / 
*                       Unaffected denominator
                        [(l1-p)^2-m^2];

L [t_gamma_3*t_gamma_f*I_Gamma] =   (g_(1,lf)+M*gi_(1))*(

(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])
+(g_(1,lf)+M*gi_(1))*(g_(1,lf)+M*gi_(1))/([lf^2-M^2]*[lf^2-M^2])

                                    )*(g_(1,lf)+M*gi_(1)) / 

                                    ([lf^2-M^2]*[lf^2-M^2]*[l2^2-m^2]*[l2^2-m^2]) /
*                                   Unaffected denominator
                                    [(l1-p)^2-m^2];

L [t_gamma_1*I_Gamma] = (g_(1,lf)+M*gi_(1))*(

(g_(1,Qf)+g_(1,l2)+M*gi_(1))*(g_(1,Qf)+g_(1,l1)+M*gi_(1))/([(Qf+l1)^2-M^2]*[(Qf+l2)^2-M^2])
+(g_(1,Qf)+g_(1,l1)-g_(1,l2)+M*gi_(1))*(g_(1,Qf)+g_(1,l1)+M*gi_(1))/([(Qf+l1-l2)^2-M^2]*[(Qf+l1)^2-M^2])
+(g_(1,Qf)+g_(1,l1)-g_(1,l2)+M*gi_(1))*(g_(1,Qf)+g_(1,Q)-g_(1,l2)+M*gi_(1))/([(Qf+l1-l2)^2-M^2]*[(Qf+Q-l2)^2-M^2])
+(g_(1,Qf)+g_(1,l2)+M*gi_(1))*(g_(1,Qf)+g_(1,Q)-g_(1,l1)+g_(1,l2)+M*gi_(1))/([(Qf+l2)^2-M^2]*[(Qf+Q-l1+l2)^2-M^2])
+(g_(1,Qf)+g_(1,Q)-g_(1,l1)+M*gi_(1))*(g_(1,Qf)+g_(1,Q)-g_(1,l2)+M*gi_(1))/([(Qf+Q-l1)^2-M^2]*[(Qf+Q-l2)^2-M^2])
+(g_(1,Qf)+g_(1,Q)-g_(1,l1)+M*gi_(1))*(g_(1,Qf)+g_(1,Q)-g_(1,l1)+g_(1,l2)+M*gi_(1))/([(Qf+Q-l1)^2-M^2]*[(Qf+Q-l1+l2)^2-M^2])
                        
                        )*(g_(1,lf)+g_(1,p)+M*gi_(1)) / 

                        ([lf^2-M^2]*[(lf+p)^2-M^2]*[(l2-l1)^2-m^2]*[(l1-Q)^2-m^2]*[l2^2-m^2]);

* Defining the t_Gamma operator.
L [t_Gamma] = (1+(p(i))*[d/dp](p)+1/2*(p(i))*[d/dp](p)*(p(j))*[d/dp](p));

Trace4,1;
.sort

#call ReplaceMomenta
.sort
hide;

* Defining the non-proper forests.
L [t_Gamma*I_Gamma] =   [t_Gamma]*[I_Gamma];
L [t_Gamma*t_gamma_f*I_Gamma] = [t_Gamma]*[t_gamma_f*I_Gamma];
L [t_Gamma*t_gamma_3*I_Gamma] = [t_Gamma]*[t_gamma_3*I_Gamma];
L [t_Gamma*t_gamma_3*t_gamma_f*I_Gamma] =   [t_Gamma]*[t_gamma_3*t_gamma_f*I_Gamma];
L [t_Gamma*t_gamma_1*I_Gamma] = [t_Gamma]*[t_gamma_1*I_Gamma];
.sort

#call Derivative

* Setting the renormalization scale to QQ=0.
id p = 0;
argument;
    id p = 0;
    id [(lf+p)^2-M^2]^n?        = [lf^2-M^2]^n;
    id [(l1-p)^2-m^2]^n?        = [l1^2-m^2]^n;
    id [(lf+p-l1)^2-M^2]^n?     = [(lf-l1)^2-M^2]^n;
    id [(lf+p-l2)^2-M^2]^n?     = [(lf-l2)^2-M^2]^n;
    id [(lf+p-l1+l2)^2-M^2]^n?  = [(lf-l1+l2)^2-M^2]^n;
endargument;
.sort
unhide;

#call RestoreMomenta
.sort

* Defining the common denominator of R_Gamma.
L [Den(R_Gamma)] = 
[lf^2-M^2]*
[(lf+p)^2-M^2]*
[(l2-l1)^2-m^2]*
[(l1-p)^2-m^2]*
[l2^2-m^2]*
[(l1-Q)^2-m^2]*
[l1^2-m^2]*

[(lf+l1)^2-M^2]*
[(lf+l2)^2-M^2]*
[(lf+p-l1)^2-M^2]*
[(lf+p-l2)^2-M^2]*
[(lf+p-l1+l2)^2-M^2]*
[(lf+l1-l2)^2-M^2]*
[(lf-l2)^2-M^2]*
[(Qf+l1)^2-M^2]*
[(Qf+l2)^2-M^2]*
[(Qf+Q-l1)^2-M^2]*
[(Qf+Q-l2)^2-M^2]*
[(Qf+Q-l1+l2)^2-M^2]*
[(Qf+l1-l2)^2-M^2]*
[(lf-l1)^2-M^2]*
[(lf-l1+l2)^2-M^2]
1
    ;

.sort
hide;

* Defining the renormalization parts over a common denominator, but NOT added together.
* Proper forests
G [I_G]            = [I_Gamma]/[Den(R_Gamma)];
G [t_g_f]          = [t_gamma_f*I_Gamma]/[Den(R_Gamma)];
G [t_g_3]          = [t_gamma_3*I_Gamma]/[Den(R_Gamma)];
G [t_g_3t_g_f]     = [t_gamma_3*t_gamma_f*I_Gamma]/[Den(R_Gamma)];
G [t_g_1]          = [t_gamma_1*I_Gamma]/[Den(R_Gamma)];

* Non-proper forests
G [t_G]            = [t_Gamma*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_f]       = [t_Gamma*t_gamma_f*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_3]       = [t_Gamma*t_gamma_3*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_3t_g_f]  = [t_Gamma*t_gamma_3*t_gamma_f*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_1]       = [t_Gamma*t_gamma_1*I_Gamma]/[Den(R_Gamma)];
.sort

* Putting everything over a common denominator
*---------------------------------------------------------------------------------------------
* Denominators of I_Gamma
id [lf^2-M^2]^k?neg_            	 = [lf^2-M^2]^(5+k)             *g([lf^2-M^2]^-4);
id [(lf+p)^2-M^2]^k?neg_        	 = [(lf+p)^2-M^2]^(2+k)         *g([(lf+p)^2-M^2]^-1);
id [(l2-l1)^2-m^2]^k?neg_       	 = [(l2-l1)^2-m^2]^(2+k)        *g([(l2-l1)^2-m^2]^-1);
id [(l1-p)^2-m^2]^k?neg_        	 = [(l1-p)^2-m^2]^(2+k)         *g([(l1-p)^2-m^2]^-1);
id [l2^2-m^2]^k?neg_            	 = [l2^2-m^2]^(3+k)             *g([l2^2-m^2]^-2);
id [(l1-Q)^2-m^2]^k?neg_        	 = [(l1-Q)^2-m^2]^(2+k)         *g([(l1-Q)^2-m^2]^-1);
id [l1^2-m^2]^k?neg_        	     = [l1^2-m^2]^(4+k)             *g([l1^2-m^2]^-3);


id [(lf+l1)^2-M^2]^k?neg_       	 = [(lf+l1)^2-M^2]^(2+k)        *g([(lf+l1)^2-M^2]^-1);
id [(lf+l2)^2-M^2]^k?neg_       	 = [(lf+l2)^2-M^2]^(3+k)        *g([(lf+l2)^2-M^2]^-2);
id [(lf+p-l1)^2-M^2]^k?neg_     	 = [(lf+p-l1)^2-M^2]^(2+k)      *g([(lf+p-l1)^2-M^2]^-1);
id [(lf+p-l2)^2-M^2]^k?neg_     	 = [(lf+p-l2)^2-M^2]^(2+k)      *g([(lf+p-l2)^2-M^2]^-1);
id [(lf+p-l1+l2)^2-M^2]^k?neg_  	 = [(lf+p-l1+l2)^2-M^2]^(2+k)   *g([(lf+p-l1+l2)^2-M^2]^-1);
id [(lf+l1-l2)^2-M^2]^k?neg_    	 = [(lf+l1-l2)^2-M^2]^(2+k)     *g([(lf+l1-l2)^2-M^2]^-1);

id [(lf-l2)^2-M^2]^k?neg_   	     = [(lf-l2)^2-M^2]^(4+k)        *g([(lf-l2)^2-M^2]^-3);
id [(Qf+l1)^2-M^2]^k?neg_       	 = [(Qf+l1)^2-M^2]^(2+k)        *g([(Qf+l1)^2-M^2]^-1);
id [(Qf+l2)^2-M^2]^k?neg_       	 = [(Qf+l2)^2-M^2]^(2+k)        *g([(Qf+l2)^2-M^2]^-1);

id [(Qf+Q-l1)^2-M^2]^k?neg_     	 = [(Qf+Q-l1)^2-M^2]^(2+k)      *g([(Qf+Q-l1)^2-M^2]^-1);
id [(Qf+Q-l2)^2-M^2]^k?neg_     	 = [(Qf+Q-l2)^2-M^2]^(2+k)      *g([(Qf+Q-l2)^2-M^2]^-1);
id [(Qf+Q-l1+l2)^2-M^2]^k?neg_  	 = [(Qf+Q-l1+l2)^2-M^2]^(2+k)   *g([(Qf+Q-l1+l2)^2-M^2]^-1);
id [(Qf+l1-l2)^2-M^2]^k?neg_    	 = [(Qf+l1-l2)^2-M^2]^(2+k)     *g([(Qf+l1-l2)^2-M^2]^-1);
id [(l2-l1)^2-m^2]^k?neg_       	 = [(l2-l1)^2-m^2]^(2+k)        *g([(l2-l1)^2-m^2]^-1);
id [(lf-l1)^2-M^2]^k?neg_     	     = [(lf-l1)^2-M^2]^(4+k)        *g([(lf-l1)^2-M^2]^-3);
id [(lf-l1+l2)^2-M^2]^k?neg_  	     = [(lf-l1+l2)^2-M^2]^(4+k)     *g([(lf-l1+l2)^2-M^2]^-3);

.sort

* Defining a procedure that writes the current progress to the terminal.
* This can then be downloaded to check the progress when running this code on a cluster.
* The procedure requires a preprocessor variable 'step' to be defined beforehand.
#define step "1"
#procedure ProgressReport
    #message Step 'step'/22
    #redefine step "{'step'+1}"
#endprocedure

*** Writing out the numerator to obtain the cancellations.
* Note that we sort after each substitution to reduce the number of terms in the intermediate expression.
#call ProgressReport
id [lf^2-M^2]                	 = (lf.lf-M^2);                                   .sort
#call ProgressReport
id [(lf+p)^2-M^2]            	 = (lf.lf+p.p+2*lf.p-M^2);                        .sort
#call ProgressReport
id [(lf+l1)^2-M^2]           	 = (lf.lf+l1.l1+2*lf.l1-M^2);                     .sort
#call ProgressReport
id [(lf+l2)^2-M^2]           	 = (lf.lf+l2.l2+2*lf.l2-M^2);                     .sort
#call ProgressReport
id [(lf+p-l1)^2-M^2]         	 = (lf.lf+p.p+l1.l1+2*lf.p-2*lf.l1-2*p.l1-M^2);   .sort
#call ProgressReport
id [(lf+p-l2)^2-M^2]         	 = (lf.lf+p.p+l2.l2+2*lf.p-2*lf.l2-2*p.l2-M^2);   .sort
#call ProgressReport
id [(lf+p-l1+l2)^2-M^2]      	 = (lf.lf+p.p+l1.l1+l2.l2
                                   +2*lf.p-2*lf.l1+2*lf.l2
                                   -2*p.l1+2*p.l2-2*l1.l2-M^2);                   .sort
#call ProgressReport
id [(lf+l1-l2)^2-M^2]        	 = (lf.lf+l1.l1+l2.l2+2*lf.l1-2*lf.l2-2*l1.l2-M^2);.sort

#call ProgressReport
id [l2^2-m^2]                	 = (l2.l2-m^2);                                   .sort
#call ProgressReport
id [(l2-l1)^2-m^2]           	 = (l1.l1+l2.l2-2*l1.l2-m^2);                     .sort
#call ProgressReport
id [(l1-p)^2-m^2]            	 = (l1.l1+p.p-2*l1.p-m^2);                        .sort

* Denominators of t_gamma_3*I_Gamma
#call ProgressReport
id [(lf-l2)^2-M^2]      	     = (lf.lf+l2.l2-2*lf.l2-M^2);
.sort

* Denominators of t_gamma_1*I_Gamma
#call ProgressReport
id [(Qf+l1)^2-M^2]           	 = (Qf.Qf+l1.l1+2*Qf.l1-M^2);
.sort
#call ProgressReport
id [(Qf+l2)^2-M^2]           	 = (Qf.Qf+l2.l2+2*Qf.l2-M^2);
.sort
#call ProgressReport
id [(Qf+Q-l1)^2-M^2]         	 = (Qf.Qf+Q.Q+l1.l1+2*Qf.Q-2*Qf.l1-2*Q.l1-M^2);
.sort
#call ProgressReport
id [(Qf+Q-l2)^2-M^2]         	 = (Qf.Qf+Q.Q+l2.l2+2*Qf.Q-2*Qf.l2-2*Q.l2-M^2);
.sort
#call ProgressReport
id [(Qf+Q-l1+l2)^2-M^2]      	 = (Qf.Qf+Q.Q+l1.l1+l2.l2
                                   +2*Qf.Q-2*Qf.l1+2*Qf.l2
                                   -2*Q.l1+2*Q.l2-2*l1.l2-M^2);
.sort
#call ProgressReport
id [(Qf+l1-l2)^2-M^2]        	 = (Qf.Qf+l1.l1+l2.l2+2*Qf.l1-2*Qf.l2-2*l1.l2-M^2);
.sort
#call ProgressReport
id [(l1-Q)^2-m^2]            	 = (l1.l1+Q.Q-2*l1.Q-m^2);  
.sort
#call ProgressReport                            

* Denominators of t_Gamma*\overline{R}_Gamma
id [l1^2-m^2]            	     = (l1.l1-m^2);
.sort
#call ProgressReport
id [(lf-l1)^2-M^2]         	     = (lf.lf+l1.l1-2*lf.l1-M^2);
.sort
#call ProgressReport
id [(lf-l1+l2)^2-M^2]      	     = (lf.lf+l1.l1+l2.l2
                                   -2*lf.l1+2*lf.l2
                                   -2*l1.l2-M^2);
.sort
#message Expansion of numerator completed!

* Restoring the denominators (removing the commuting function 'g').
id g(a?) = a;
.store
#message Denominators restored!

#message Saving results...
save Beetle_Diagram_Forests.sav;
Delete storage;
.end