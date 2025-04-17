% make G with regularized stokeslet
function G=makeGregstv2(sp,mu,epsvec)    % epsvec : 1 row x N column of e,  
nsp=length(sp) ;                         % sp: 3 rows x N columns of position matrix (r)
G = zeros(3*nsp) ;                       % 3 = 3D, nsp = length of r
L1 = zeros(1,3*nsp) ;                    % G matrix describes each sphere in x,y,z position in relation to all other spheres
L2 = L1 ; 
L3 = L1 ; 
in1 = 1:3:3*nsp-2 ; 
in2 = 2:3:3*nsp-1 ; 
in3 = 3:3:3*nsp ; 
term1 = zeros(3,3*nsp) ; 
coefs = ones(3,3*nsp) ; 


% seqterm1 = 1:9:9*nsp-8 ; 
% seqtot_term1 = zeros(1,3*nsp) ; 
% seqtot_term1 = [seqterm1,seqterm1+4,seqterm1+8] ; 
% seqtot_term1(1,in2) = seqterm1+4 ; 
% seqtot_term1(1,in3) = seqterm1+8 ; 



for i=1:nsp
    r_ab = repmat(sp(:,i),1,nsp)-sp ; 
    nr_ab = vecnorm(r_ab) ; 
    nrab2 = nr_ab.^2;
    nrab2eps2 = nrab2+epsvec.^2 ;
    nrab22eps2 = nrab2+2*epsvec.^2 ;
    L1(in1) = (r_ab(1,:).^2 + nrab22eps2) ; 
    L1(in2) = (r_ab(1,:).*r_ab(2,:))  ; 
    L1(in3) = (r_ab(1,:).*r_ab(3,:))  ; 
    L2(in1) = (L1(in2))  ; 
    L2(in2) = (r_ab(2,:).^2 + nrab22eps2)   ;
    L2(in3) = (r_ab(2,:).*r_ab(3,:))  ;
    L3(in1) = (L1(in3))  ; 
    L3(in2) = (L2(in3))  ; 
    L3(in3) = (r_ab(3,:).^2 + nrab22eps2)   ;
    L123 = [L1;L2;L3] ; 
     
%     term1(1,in1) = nrab2eps2 ; 
%     term1(2,in2) = nrab2eps2 ; 
%     term1(3,in3) = nrab2eps2 ; 
    
%     term1(seqtot_term1) = repmat(nrab2eps2,1,3) ; 
    
    coefs(1,in1) = nrab2eps2 ; 
    coefs(1,in2) = nrab2eps2 ; 
    coefs(1,in3) = nrab2eps2 ; 
    coefs(2,:) = coefs(1,:) ; 
    coefs(3,:) = coefs(1,:) ; 
    
%     size(nrab2eps2)
%     size(L123)
%     size(term1)
%     G(3*i-2:3*i,:) = (1/(mu*8*pi)./coefs.^1.5).*( term1 + L123 ) ; 
    G(3*i-2:3*i,:) = (1/(mu*8*pi)./coefs.^1.5).*L123 ; 

end
