;  procedure sinfunc
;
;Procedure to produce a function A(0)*sin^(A(1))(x)
;
;inputs:   X:  the value (angle) X in degrees
;	   A:  the input parameters for the function
;outputs;  F:  the result
;
PRO sinfunc, X, A, F

	F = A(0) * (sin(X*!pi/180))^A(1)

END