%-------------------------------------------------------------------------%
function [xxx] = VecElementInterp (xx1,xx2, min, sec)
%-------------------------------------------------------------------------%

min = min/60;
sec = sec/60/60;

y1 = [xx1(1); xx2(1)];
y2 = [xx1(2); xx2(2)];
y3 = [xx1(3); xx2(3)];
x = [0; 1];
xv = min+sec;

interpResult1 = interp1(x,y1,xv);
interpResult2 = interp1(x,y2,xv);
interpResult3 = interp1(x,y3,xv);

xxx = [interpResult1; interpResult2; interpResult3];
