* This FORM code performs a UV power counting of the Rocket Diagram in the different regions of integration space.
* This code loads the expressions for the Forests of the Rocket from the same directory as where this code is located.
* These expressions can be generated using the file "Rocket_Diagram_Forests.frm".
* Note that the contribution from the measure has NOT been taken into account and should be corrected for manually.

Symbols m,M;
Symbols a,n;
Symbols [l1^2-M^2],[(l1-l)^2-M^2],[(l1+p1+p2)^2-M^2],[(l1+p2)^2-M^2],
        [(l1+l)^2-M^2],[(l1+p1+l)^2-M^2],[(l1-p2)^2-M^2],
        [(l+p1+p2)^2-m^2],[l^2-m^2],
        [l2^2-M^2],[(l2+p1+p2)^2-M^2],[(l2-l)^2-M^2];       * Denominators of the Rocket
Vectors p1,p2,l,l1,l2;      * unfixed momenta
CFunctions F (Symmetric);

* Turning warnings off to prevent FORM from giving one every time a term does not fit inside the bracket function "F" (which happens all the time).
OFF Warnings;
OFF Statistics;

.global

* Load the previously computed expressions for the BPHZ renormalized Rocket diagram forests.
* ---------------------------------------------------------------------------------
load Rocket_Diagram_Forests.sav;
* ---------------------------------------------------------------------------------

* Performing the whole power counting procedure for a variety of forests and sums thereof.
* This allows us to track how the behaviour changes with subtractions and deduce general rules from those changes.
#do B = {[I_G]-[t_g_1]-[t_g_2]+[t_g_1t_g_2]-[t_G]+[t_Gt_g_1]+[t_Gt_g_2]-[t_Gt_g_1t_g_2]}
#message Power counting for `B'


G [dummy] =
`B'
;

.sort
hide [dummy];

* Scaling l by a factor 'a' to do the power counting in the region {l}.
L [R_Gamma] = [dummy];      * A quick fix. Without these two lines, all power countings in {l} will return 0. No idea why.

* Removing all odd powers of loop momenta
* ---------------------------------------
* #do C = {l,l1,l2}
* #message Removing odd powers of `C'.
*     id `C' = a*`C';
*     id a^n?odd_ = 0;
*     id a = 1;
*     .sort
* #enddo
* ---------------------------------------

.sort
hide;
L [count l] = [R_Gamma];
id l = a*l;
id [(l1-l)^2-M^2]^n?        = (a^2)^n;
id [(l1+l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+l)^2-M^2]^n?     = (a^2)^n;
id [(l+p1+p2)^2-m^2]^n?     = (a^2)^n;
id [l^2-m^2]^n?             = (a^2)^n;
id [(l2-l)^2-M^2]^n?        = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l];



* Scaling l1 by a factor 'a' to do the power counting in the region {l1}.
L [count l1] = [R_Gamma];
id l1 = a*l1;
id [l1^2-M^2]^n?            = (a^2)^n;
id [(l1-l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l1+p2)^2-M^2]^n?       = (a^2)^n;
id [(l1+l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+l)^2-M^2]^n?     = (a^2)^n;
id [(l1-p2)^2-M^2]^n?       = (a^2)^n;

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
id [l2^2-M^2]^n?            = (a^2)^n;
id [(l2+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l2-l)^2-M^2]^n?        = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l2];



* Scaling l and l1 by a factor 'a' to do the power counting in the region {l,l1}.
L [count ll1] = [R_Gamma];
id l = a*l;
id l1 = a*l1;
id [l1^2-M^2]^n?            = (a^2)^n;
id [(l1-l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l1+p2)^2-M^2]^n?       = (a^2)^n;
id [(l1+l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+l)^2-M^2]^n?     = (a^2)^n;
id [(l1-p2)^2-M^2]^n?       = (a^2)^n;
id [(l+p1+p2)^2-m^2]^n?     = (a^2)^n;
id [l^2-m^2]^n?             = (a^2)^n;
id [(l2-l)^2-M^2]^n?        = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count ll1];


* Scaling l and l2 by a factor 'a' to do the power counting in the region {l,l2}.
L [count ll2] = [R_Gamma];
id l = a*l;
id l2 = a*l2;
id [(l1-l)^2-M^2]^n?        = (a^2)^n;
id [(l1+l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+l)^2-M^2]^n?     = (a^2)^n;
id [(l+p1+p2)^2-m^2]^n?     = (a^2)^n;
id [l^2-m^2]^n?             = (a^2)^n;
id [l2^2-M^2]^n?            = (a^2)^n;
id [(l2+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l2-l)^2-M^2]^n?        = (a^2)^n;

* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count ll2];


* Scaling l1 and l2 by a factor 'a' to do the power counting in the region {l1,l2}.
L [count l1l2] = [R_Gamma];
id l1 = a*l1;
id l2 = a*l2;
id [l1^2-M^2]^n?            = (a^2)^n;
id [(l1-l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l1+p2)^2-M^2]^n?       = (a^2)^n;
id [(l1+l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+l)^2-M^2]^n?     = (a^2)^n;
id [(l1-p2)^2-M^2]^n?       = (a^2)^n;
id [l2^2-M^2]^n?            = (a^2)^n;
id [(l2+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l2-l)^2-M^2]^n?        = (a^2)^n;


* Factoring out the powers of 'a' to read off the power counting.
Bracket a;
.sort
* Put all remaining polynomial terms inside a function F and set this to 1 (effectively throwing them away)
Collect F;
id F(a?) = 1;

Print;
.sort
hide [count l1l2];


* Scaling l1, l2 and l by a factor 'a' to do the power counting in the region {l,l1,l2}.
L [count ll1l2] = [R_Gamma];
id l = a*l;
id l1 = a*l1;
id l2 = a*l2;
id [l1^2-M^2]^n?            = (a^2)^n;
id [(l1-l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l1+p2)^2-M^2]^n?       = (a^2)^n;
id [(l1+l)^2-M^2]^n?        = (a^2)^n;
id [(l1+p1+l)^2-M^2]^n?     = (a^2)^n;
id [(l1-p2)^2-M^2]^n?       = (a^2)^n;
id [(l+p1+p2)^2-m^2]^n?     = (a^2)^n;
id [l^2-m^2]^n?             = (a^2)^n;
id [l2^2-M^2]^n?            = (a^2)^n;
id [(l2+p1+p2)^2-M^2]^n?    = (a^2)^n;
id [(l2-l)^2-M^2]^n?        = (a^2)^n;

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