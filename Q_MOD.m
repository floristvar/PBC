function [Qb,Sb] = Q_MOD(theta,S1)
% theta: the angle of the inclination of composite lamina
sa=sind(theta); 
ca=cosd(theta);
% estimating compliance matrix elemnents in 1-2 coordinates
sb11 = S1(1,1)*ca^4 + (2*S1(1,2)+S1(3,3))*sa^2*ca^2 + S1(2,2)*sa^4;
sb12 = S1(1,2)*(sa^4+ca^4) + (S1(1,1)+S1(2,2)-S1(3,3))*sa^2*ca^2;
sb22 = S1(1,1)*sa^4 + (2*S1(1,2)+S1(3,3))*sa^2*ca^2 + S1(2,2)*ca^4;
sb16 = (2*S1(1,1)-2*S1(1,2)-S1(3,3))*sa*ca^3 - (2*S1(2,2)-2*S1(1,2)-S1(3,3))*sa^3*ca;
sb26 = (2*S1(1,1)-2*S1(1,2)-S1(3,3))*sa^3*ca - (2*S1(2,2)-2*S1(1,2)-S1(3,3))*sa*ca^3;
sb66 = 2*(2*S1(1,1)+2*S1(2,2)-4*S1(1,2)-S1(3,3))*sa^2*ca^2 + S1(3,3)*(sa^4+ca^4);
%final transfomed matrix Sbar matrix 
Sb =[sb11 sb12 sb16; sb12 sb22 sb26; sb16 sb26 sb66];
Qb=inv(Sb);
end