* This code computes the expressions for each term in the Forest Formula for the Kite Diagram and puts these over a common denominator so they can be added together later.
* It stores the resulting expressions in a file "Kite_Diagram_Forests.sav", which is placed in the same directory as this code.

Dimension 4;
Symbols a,n,m,M;
Symbols [l1^2-M^2],[(l1+p)^2-M^2],[(l1-ls)^2-M^2],
        [l2^2-M^2],[(l2+p)^2-M^2],[(l2-ls)^2-M^2],
        [ls^2-m^2],[(ls+p)^2-m^2];                              * Denominators
Functions [d/dp],f (Symmetric);
CFunctions g (Symmetric);
Vectors x,y,z;
Vectors p,ls,l1,l2;         * Unfixed momenta
Indices i,j;

OFF Statistics;

*** Definitions of procedures we will need later.
#procedure FactorDenominators
Bracket 
    [l1^2-M^2],[(l1+p)^2-M^2],[(l1-ls)^2-M^2],
    [l2^2-M^2],[(l2+p)^2-M^2],[(l2-ls)^2-M^2],
    [ls^2-m^2],[(ls+p)^2-m^2];
#endprocedure

#define step "1"
* Defining a procedure that writes the current progress to the terminal.
* This can then be downloaded to check the progress when running this code on a cluster.
* The procedure requires a preprocessor variable 'step' to be defined beforehand.
* To use it: call the procedure BEFORE each noteworthy (computationally expensive) step.
#procedure ProgressReport
    #message Step 'step'/8
    #redefine step "{'step'+1}"
#endprocedure

#procedure Derivative
Repeat;
* All occurrences of dot products.
    id [d/dp](x?)*f(y?!{x?},z?!{x?})       = f(y,z)*[d/dp](x);
    id [d/dp](x?)*f(z?!{x?},g(y?))         = f(g(y),z)*[d/dp](x);
    id [d/dp](x?)*f(g(y?),g(z?))           = f(g(y),g(z))*[d/dp](x);

    id g(x?(i?))*[d/dp](x?)*f(x?,x?)                 = 2*f(g(x),x)+f(x,x)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(y?!{x?},x?)            = g(g(x),y)+f(y,x)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(a?,x?)                 = g(g(x),a)+f(a,x)*g(x(i))*[d/dp](x);
    id g(x?(i?))*[d/dp](x?)*f(g(y?),x?)              = g(g(x),g(y))+f(g(y),x)*g(x(i))*[d/dp](x);

* All denominators.
    id g(x?{l1}(i?))*[d/dp](x?{l1})*f([l1^2-M^2]^n?)          = (+1*n*2*(f(l1,g(x)))*f([l1^2-M^2]^(n-1)) + f([l1^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{l1}(i?))*[d/dp](x?!{l1})*f([l1^2-M^2]^n?)        = (f([l1^2-M^2]^n)*g(x(i))*[d/dp](x));

    id g(x?{l1,p}(i?))*[d/dp](x?{l1,p})*f([(l1+p)^2-M^2]^n?)  = (+1*n*2*(f(l1,g(x))+f(p,g(x)))*f([(l1+p)^2-M^2]^(n-1)) + f([(l1+p)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{l1,p}(i?))*[d/dp](x?!{l1,p})*f([(l1+p)^2-M^2]^n?)= (f([(l1+p)^2-M^2]^n)*g(x(i))*[d/dp](x));

    id g(x?{l1}(i?))*[d/dp](x?{l1})*f([(l1-ls)^2-M^2]^n?)     = (+1*n*2*(f(l1,g(x))-f(ls,g(x)))*f([(l1-ls)^2-M^2]^(n-1)) + f([(l1-ls)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?{ls}(i?))*[d/dp](x?{ls})*f([(l1-ls)^2-M^2]^n?)     = (-1*n*2*(f(l1,g(x))-f(ls,g(x)))*f([(l1-ls)^2-M^2]^(n-1)) + f([(l1-ls)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{l1,ls}(i?))*[d/dp](x?!{l1,ls})*f([(l1-ls)^2-M^2]^n?)   = (f([(l1-ls)^2-M^2]^n)*g(x(i))*[d/dp](x));


    id g(x?{l2}(i?))*[d/dp](x?{l2})*f([l2^2-M^2]^n?)          = (+1*n*2*(f(l2,g(x)))*f([l2^2-M^2]^(n-1)) + f([l2^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{l2}(i?))*[d/dp](x?!{l2})*f([l2^2-M^2]^n?)        = (f([l2^2-M^2]^n)*g(x(i))*[d/dp](x));

    id g(x?{l2,p}(i?))*[d/dp](x?{l2,p})*f([(l2+p)^2-M^2]^n?)  = (+1*n*2*(f(l2,g(x))+f(p,g(x)))*f([(l2+p)^2-M^2]^(n-1)) + f([(l2+p)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{l2,p}(i?))*[d/dp](x?!{l2,p})*f([(l2+p)^2-M^2]^n?)= (f([(l2+p)^2-M^2]^n)*g(x(i))*[d/dp](x));

    id g(x?{l2}(i?))*[d/dp](x?{l2})*f([(l2-ls)^2-M^2]^n?)     = (+1*n*2*(f(l2,g(x))-f(ls,g(x)))*f([(l2-ls)^2-M^2]^(n-1)) + f([(l2-ls)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?{ls}(i?))*[d/dp](x?{ls})*f([(l2-ls)^2-M^2]^n?)     = (-1*n*2*(f(l2,g(x))-f(ls,g(x)))*f([(l2-ls)^2-M^2]^(n-1)) + f([(l2-ls)^2-M^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{l2,ls}(i?))*[d/dp](x?!{l2,ls})*f([(l2-ls)^2-M^2]^n?)   = (f([(l2-ls)^2-M^2]^n)*g(x(i))*[d/dp](x));


    id g(x?{ls}(i?))*[d/dp](x?{ls})*f([ls^2-m^2]^n?)          = (+1*n*2*(f(ls,g(x)))*f([ls^2-m^2]^(n-1)) + f([ls^2-m^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{ls}(i?))*[d/dp](x?!{ls})*f([ls^2-m^2]^n?)        = (f([ls^2-m^2]^n)*g(x(i))*[d/dp](x));

    id g(x?{ls,p}(i?))*[d/dp](x?{ls,p})*f([(ls+p)^2-m^2]^n?)  = (+1*n*2*(f(ls,g(x))+f(p,g(x)))*f([(ls+p)^2-m^2]^(n-1)) + f([(ls+p)^2-m^2]^n)*g(x(i))*[d/dp](x));
    id g(x?!{ls,p}(i?))*[d/dp](x?!{ls,p})*f([(ls+p)^2-m^2]^n?)= (f([(ls+p)^2-m^2]^n)*g(x(i))*[d/dp](x));
EndRepeat;
id [d/dp](x?) = 0;
#endprocedure

#procedure ReplaceMomenta
* Replacing mometa by noncommuting functions so they don't commute with the derivative operator.
    id x?.y? = f(x,y);
    id [l1^2-M^2]^n?            = f([l1^2-M^2]^n);
    id [(l1+p)^2-M^2]^n?        = f([(l1+p)^2-M^2]^n);
    id [(l1-ls)^2-M^2]^n?       = f([(l1-ls)^2-M^2]^n);
    id [l2^2-M^2]^n?            = f([l2^2-M^2]^n);
    id [(l2+p)^2-M^2]^n?        = f([(l2+p)^2-M^2]^n);
    id [(l2-ls)^2-M^2]^n?       = f([(l2-ls)^2-M^2]^n);
    id [ls^2-m^2]^n?            = f([ls^2-m^2]^n);
    id [(ls+p)^2-m^2]^n?        = f([(ls+p)^2-m^2]^n);
    id f(1) = 1;
    id x?(i?) = g(x(i));
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

*** Defining the forests of the Kite.
L [t_Gamma] = (1+(p(i))*[d/dp](p)+1/2*(p(i))*[d/dp](p)*(p(j))*[d/dp](p));

* Atomic forests.
L [I_Gamma] =   ((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,ls)+M*gi_(1))*(g_(1,l1)+g_(1,p)+M*gi_(1))) / 
                ([l1^2-M^2]*[(l1+p)^2-M^2]*[(l1-ls)^2-M^2])*

                1/([ls^2-m^2]*[(ls+p)^2-m^2])*

                ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+g_(2,p)+M*gi_(2))*(g_(2,l2)-g_(2,ls)+M*gi_(2))) / 
                ([l2^2-M^2]*[(l2+p)^2-M^2]*[(l2-ls)^2-M^2]);

L [t_gamma_1*I_Gamma] = ((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                        ([l1^2-M^2]^3)*

                        1/([ls^2-m^2]*[(ls+p)^2-m^2])*

                        ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+g_(2,p)+M*gi_(2))*(g_(2,l2)-g_(2,ls)+M*gi_(2))) / 
                        ([l2^2-M^2]*[(l2+p)^2-M^2]*[(l2-ls)^2-M^2]);

L [t_gamma_2*I_Gamma] = ((g_(1,l1)+M*gi_(1))*(g_(1,l1)-g_(1,ls)+M*gi_(1))*(g_(1,l1)+g_(1,p)+M*gi_(1))) / 
                        ([l1^2-M^2]*[(l1+p)^2-M^2]*[(l1-ls)^2-M^2])*

                        1/([ls^2-m^2]*[(ls+p)^2-m^2])*

                        ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))) / 
                        ([l2^2-M^2]^3);
      
* Composite forests
L [t_gamma_1*t_gamma_2*I_Gamma] =   ((g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))*(g_(1,l1)+M*gi_(1))) / 
                                    ([l1^2-M^2]^3)*

                                    1/([ls^2-m^2]*[(ls+p)^2-m^2])*

                                    ((g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))*(g_(2,l2)+M*gi_(2))) / 
                                    ([l2^2-M^2]^3);

Trace4,1;
Trace4,2;
.sort

#call ReplaceMomenta
.sort
hide;

* Non-proper forests.
L [t_Gamma*I_Gamma] = [t_Gamma]*[I_Gamma];
L [t_Gamma*t_gamma_1*I_Gamma] = [t_Gamma]*[t_gamma_1*I_Gamma];
L [t_Gamma*t_gamma_2*I_Gamma] = [t_Gamma]*[t_gamma_2*I_Gamma];

L [t_Gamma*t_gamma_1*t_gamma_2*I_Gamma] = [t_Gamma]*[t_gamma_1*t_gamma_2*I_Gamma];

#call Derivative

* Set the renormalization scale.
id p = 0;
id f(p,x?) = f(0,x);
id f(p,g(x?)) = f(0,g(x));
argument;
    id [(l1+p)^2-M^2]^n? = [l1^2-M^2]^n;
    id [(l2+p)^2-M^2]^n? = [l2^2-M^2]^n;
    id [(ls+p)^2-m^2]^n? = [ls^2-m^2]^n;
endargument;
.sort
unhide;

#call RestoreMomenta
.sort

* Defining the common denominator of R_Gamma.
L [Den(R_Gamma)] = 
    [l1^2-M^2]*[(l1+p)^2-M^2]*[(l1-ls)^2-M^2]*
    [l2^2-M^2]*[(l2+p)^2-M^2]*[(l2-ls)^2-M^2]*
    [ls^2-m^2]*[(ls+p)^2-m^2]
;

.sort
hide;

* Defining the forests that are put over a common denominator (abbreviated notation so FORM won't complain about the names being too long just like this comment).
G [I_G] = [I_Gamma]/[Den(R_Gamma)];
G [t_g_1] = [t_gamma_1*I_Gamma]/[Den(R_Gamma)];
G [t_g_2] = [t_gamma_2*I_Gamma]/[Den(R_Gamma)];

G [t_g_1t_g_2] = [t_gamma_1*t_gamma_2*I_Gamma]/[Den(R_Gamma)];


G [t_G] = [t_Gamma*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_1] = [t_Gamma*t_gamma_1*I_Gamma]/[Den(R_Gamma)];
G [t_Gt_g_2] = [t_Gamma*t_gamma_2*I_Gamma]/[Den(R_Gamma)];

G [t_Gt_g_1t_g_2] = [t_Gamma*t_gamma_1*t_gamma_2*I_Gamma]/[Den(R_Gamma)];
.sort

* Putting everything over a common denominator.
id [l1^2-M^2]^n?neg_            = [l1^2-M^2]^(n+5)      * g([l1^2-M^2]^-4);
id [(l1+p)^2-M^2]^n?neg_        = [(l1+p)^2-M^2]^(n+2)  * g([(l1+p)^2-M^2]^-1);
id [(l1-ls)^2-M^2]^n?neg_       = [(l1-ls)^2-M^2]^(n+2) * g([(l1-ls)^2-M^2]^-1);
id [l2^2-M^2]^n?neg_            = [l2^2-M^2]^(n+5)      * g([l2^2-M^2]^-4);
id [(l2+p)^2-M^2]^n?neg_        = [(l2+p)^2-M^2]^(n+2)  * g([(l2+p)^2-M^2]^-1);
id [(l2-ls)^2-M^2]^n?neg_       = [(l2-ls)^2-M^2]^(n+2) * g([(l2-ls)^2-M^2]^-1);
id [ls^2-m^2]^n?neg_            = [ls^2-m^2]^(n+5)      * g([ls^2-m^2]^-4);
id [(ls+p)^2-m^2]^n?neg_        = [(ls+p)^2-m^2]^(n+2)  * g([(ls+p)^2-m^2]^-1);
.sort

* Writing out the numerator to obtain the cancellations.
* Note that we sort after each substitution to reduce the number of terms in the intermediate expression.
#call ProgressReport
id [l1^2-M^2]           = (l1.l1-M^2);
.sort
#call ProgressReport
id [(l1+p)^2-M^2]       = (l1.l1+p.p+2*l1.p-M^2);
.sort
#call ProgressReport
id [(l1-ls)^2-M^2]      = (l1.l1+ls.ls-2*l1.ls-M^2);
.sort
#call ProgressReport
id [l2^2-M^2]           = (l2.l2-M^2);
.sort
#call ProgressReport
id [(l2+p)^2-M^2]       = (l2.l2+p.p+2*l2.p-M^2);
.sort
#call ProgressReport
id [(l2-ls)^2-M^2]      = (l2.l2+ls.ls-2*l2.ls-M^2);
.sort
#call ProgressReport
id [ls^2-m^2]           = (ls.ls-m^2);
.sort
#call ProgressReport
id [(ls+p)^2-m^2]       = (ls.ls+p.p+2*ls.p-m^2);
.sort

* Restoring the denominators (removing the noncommuting function 'g').
id g(a?) = a;

#call FactorDenominators

.store
save Kite_Diagram_Forests.sav;
Delete storage;
.end
