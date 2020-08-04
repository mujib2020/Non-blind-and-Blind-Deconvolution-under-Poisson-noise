function psnr=PSNR(ud,u)
[x,y]=size(ud);
o=sum(sum( (ud-u).^2 ));
t=o/(x*y);
maxu=max(ud(:));
psnr=10*log10(maxu^2/t);
end