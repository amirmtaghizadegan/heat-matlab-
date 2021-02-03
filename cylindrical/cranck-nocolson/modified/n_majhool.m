function x=n_majhool(zarayeb,javab)

za=det(zarayeb);
x=zeros(1,length(zarayeb));
for i=1:length(zarayeb)
b=zarayeb;
b(:,i)=javab;
zb=det(b);
x(i)=zb/za;
end