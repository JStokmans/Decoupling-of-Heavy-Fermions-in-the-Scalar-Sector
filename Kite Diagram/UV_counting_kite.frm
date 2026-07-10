* This FORM code performs a UV power counting of the Kite Diagram in the different regions of integration space.
* This code loads the expressions for the Forests of the Kite from the same directory as where this code is located.
* These expressions can be generated using the file "Kite_Diagram_Forests.frm".
* Note that the contribution from the measure has NOT been taken into account and should be corrected for manually.

Dimension 4;
Symbols a,n,m,M;
Symbols [l1^2-M^2],[(l1+p)^2-M^2],[(l1-ls)^2-M^2],
        [l2^2-M^2],[(l2+p)^2-M^2],[(l2-ls)^2-M^2],
        [ls^2-m^2],[(ls+p)^2-m^2];                              * Denominators
Symbols [(ls+Q)^2-m^2],[(Qf-ls)^2-M^2],
        [(l1+Q)^2-M^2],[(l2+Q)^2-M^2];
Functions [d/dp],f (Symmetric);
CFunctions g (Symmetric);
Vectors x,y,z;
Vectors p,ls,l1,l2;         * Unfixed momenta
Vectors Q,Qf;               * phi->psi-bar+psi renormalization scale
CFunctions F (Symmetric);

* Turning warnings off to prevent FORM from giving one every time a term does not fit inside the bracket function "F" (which happens all the time).
OFF Warnings;
OFF Statistics;

.global

* Load the previously computed expressions for the BPHZ renormalized Kite diagram forests.
load Kite_Diagram_Forests.sav;

* Performing the whole power counting procedure for a variety of forests and sums thereof.
* This allows us to track how the behaviour changes with subtractions and deduce general rules from those changes.
#do B= {[I_G]-[t_g_1]-[t_g_2]+[t_g_1t_g_2]-[t_G]+[t_Gt_g_1]+[t_Gt_g_2]-[t_Gt_g_1t_g_2]}
#message Power counting for `B'

G [R_Gamma] =
`B'
;

.sort
hide [R_Gamma];

* Scaling ls by a factor 'a' to do the power counting in the region {ls}.
L [dummy] = [R_Gamma];      * A quick fix. Without these two lines, all power countings in {ls} will return 0. No idea why.
.sort
hide;
L [count ls] = [R_Gamma];
id ls = a*ls;
id [(l1-ls)^2-M^2]^n?   = (a^2)^n;
id [ls^2-m^2]^n?        = (a^2)^n;
id [(ls+p)^2-m^2]^n?    = (a^2)^n;
id [(l2-ls)^2-M^2]^n?   = (a^2)^n;
id [(ls+Q)^2-m^2]^n?    = (a^2)^n;
id [(Qf-ls)^2-M^2]^n?   = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count ls];



* Scaling l1 by a factor 'a' to do the power counting in the region {l1}.
L [count l1] = [R_Gamma];
id l1 = a*l1;
id [l1^2-M^2]^n?        = (a^2)^n;
id [(l1+p)^2-M^2]^n?    = (a^2)^n;
id [(l1-ls)^2-M^2]^n?   = (a^2)^n;
id [(l1+Q)^2-M^2]^n?    = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l1];



* Scaling l2 by a factor 'a' to do the power counting in the region {l2}.
L [count l2] = [R_Gamma];
id l2 = a*l2;
id [l2^2-M^2]^n?        = (a^2)^n;
id [(l2+p)^2-M^2]^n?    = (a^2)^n;
id [(l2-ls)^2-M^2]^n?   = (a^2)^n;
id [(l2+Q)^2-M^2]^n?    = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l2];



* Scaling ls and l1 by a factor 'a' to do the power counting in the region {ls,l1}.
L [count lsl1] = [R_Gamma];
id ls = a*ls;
id l1 = a*l1;
id [l1^2-M^2]^n?        = (a^2)^n;
id [(l1+p)^2-M^2]^n?    = (a^2)^n;
id [(l1-ls)^2-M^2]^n?   = (a^2)^n;
id [ls^2-m^2]^n?        = (a^2)^n;
id [(ls+p)^2-m^2]^n?    = (a^2)^n;
id [(l2-ls)^2-M^2]^n?   = (a^2)^n;
id [(ls+Q)^2-m^2]^n?    = (a^2)^n;
id [(Qf-ls)^2-M^2]^n?   = (a^2)^n;
id [(l1+Q)^2-M^2]^n?    = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count lsl1];


* Scaling ls and l2 by a factor 'a' to do the power counting in the region {ls,l2}.
L [count lsl2] = [R_Gamma];
id ls = a*ls;
id l2 = a*l2;
id [(l1-ls)^2-M^2]^n?   = (a^2)^n;
id [ls^2-m^2]^n?        = (a^2)^n;
id [(ls+p)^2-m^2]^n?    = (a^2)^n;
id [l2^2-M^2]^n?        = (a^2)^n;
id [(l2+p)^2-M^2]^n?    = (a^2)^n;
id [(l2-ls)^2-M^2]^n?   = (a^2)^n;
id [(ls+Q)^2-m^2]^n?    = (a^2)^n;
id [(Qf-ls)^2-M^2]^n?   = (a^2)^n;
id [(l2+Q)^2-M^2]^n?    = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count lsl2];


* Scaling l1 and l2 by a factor 'a' to do the power counting in the region {l1,l2}.
L [count l1l2] = [R_Gamma];
id l1 = a*l1;
id l2 = a*l2;
* id [l1^2-M^2]^n?        = (a^2)^n;
* id [(l1+p)^2-M^2]^n?    = (a^2)^n;
* id [(l1-ls)^2-M^2]^n?   = (a^2)^n;
* id [l2^2-M^2]^n?        = (a^2)^n;
* id [(l2+p)^2-M^2]^n?    = (a^2)^n;
* id [(l2-ls)^2-M^2]^n?   = (a^2)^n;
* id [(l1+Q)^2-M^2]^n?    = (a^2)^n;
* id [(l2+Q)^2-M^2]^n?    = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l1l2];


* Scaling l1, l2 and ls by a factor 'a' to do the power counting in the region {ls,l1,l2}.
L [count lsl1l2] = [R_Gamma];
id ls = a*ls;
id l1 = a*l1;
id l2 = a*l2;
id [l1^2-M^2]^n?        = (a^2)^n;
id [(l1+p)^2-M^2]^n?    = (a^2)^n;
id [(l1-ls)^2-M^2]^n?   = (a^2)^n;
id [ls^2-m^2]^n?        = (a^2)^n;
id [(ls+p)^2-m^2]^n?    = (a^2)^n;
id [l2^2-M^2]^n?        = (a^2)^n;
id [(l2+p)^2-M^2]^n?    = (a^2)^n;
id [(l2-ls)^2-M^2]^n?   = (a^2)^n;
id [(ls+Q)^2-m^2]^n?    = (a^2)^n;
id [(Qf-ls)^2-M^2]^n?   = (a^2)^n;
id [(l1+Q)^2-M^2]^n?    = (a^2)^n;
id [(l2+Q)^2-M^2]^n?    = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide;

#enddo
Delete storage;
.end
