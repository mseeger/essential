num=50;
maxn=700; maxm=700;

%oldnum=num; num=0; % SWITCH OFF

% DGEMM
fprintf(1,'DGEMM tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  k=1+floor(rand*maxn);
  alpha=rand; beta=rand;
  tra=double(rand>0.5); trb=double(rand>0.5);
  c=fstdebug_sampmat(m,n);
  if tra
    a=fstdebug_sampmat(k,m);
  else
    a=fstdebug_sampmat(m,k);
  end
  if trb
    b=fstdebug_sampmat(n,k);
  else
    b=fstdebug_sampmat(k,n);
  end
  u=rand;
  if u<(1/3)
    ret2=fstdeb_dgemm(c,a,tra,b,trb);
    ret1=c(:,:);
    fst_dgemm(ret1,a,tra,b,trb);
    alpha=1; beta=0;
  elseif u<(2/3)
    ret2=fstdeb_dgemm(c,a,tra,b,trb,alpha);
    ret1=c(:,:);
    fst_dgemm(ret1,a,tra,b,trb,alpha);
    beta=0;
  else
    ret2=fstdeb_dgemm(c,a,tra,b,trb,alpha,beta);
    ret1=c(:,:);
    fst_dgemm(ret1,a,tra,b,trb,alpha,beta);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dgemm_error','c','a','b','tra','trb','alpha','beta', ...
	 'ret1','ret2');
    error('ERROR (DGEMM)!');
  end
end
fprintf(1,'OK, done.\n');

% DIAGMUL
fprintf(1,'DIAGMUL tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  trs=double(rand>0.5);
  d=[];
  a=fstdebug_sampmat(m,n);
  b=fstdebug_sampmat(m,n);
  if rand>0.5
    if trs
      d=2*rand(n,1)-1;
    else
      d=2*rand(m,1)-1;
    end
  end
  if isempty(d)
    ret1=fst_diagmul(a,b,trs);
    ret2=fstdeb_diagmul(a,b,trs);
  else
    ret1=fst_diagmul(a,b,trs,d);
    ret2=fstdeb_diagmul(a,b,trs,d);
  end
  if max(max(reldiff(ret1,ret2,1e-15)))>(1e-7)
    save('diagmul_error','a','b','trs','d','ret1','ret2');
    error('ERROR (MULDIAG)!');
  end
end
fprintf(1,'OK, done.\n');

% DSYMM
fprintf(1,'DSYMM tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  left=double(rand>0.5);
  alpha=rand; beta=rand;
  c=fstdebug_sampmat(m,n);
  if left
    a=fstdebug_sampmat(m,m,1);
    b=fstdebug_sampmat(m,n);
    if rand>0.5
      a{3}='L ';
    else
      a{3}='U ';
    end
  else
    a=fstdebug_sampmat(m,n);
    b=fstdebug_sampmat(n,n,1);
    if rand>0.5
      b{3}='L ';
    else
      b{3}='U ';
    end
  end
  u=rand;
  if u<(1/3)
    ret2=fstdeb_dsymm(c,a,b);
    ret1=c(:,:);
    fst_dsymm(ret1,a,b);
    alpha=1; beta=0;
  elseif u<(2/3)
    ret2=fstdeb_dsymm(c,a,b,alpha);
    ret1=c(:,:);
    fst_dsymm(ret1,a,b,alpha);
    beta=0;
  else
    ret2=fstdeb_dsymm(c,a,b,alpha,beta);
    ret1=c(:,:);
    fst_dsymm(ret1,a,b,alpha,beta);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dsymm_error','c','a','b','alpha','beta','ret1', ...
	 'ret2');
    error('ERROR (DSYMM)!');
  end
end
fprintf(1,'OK, done.\n');

% DSYRK
fprintf(1,'DSYRK tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  trs=double(rand>0.5);
  alpha=rand; beta=rand;
  c=fstdebug_sampmat(m,m,1);
  if rand>0.5
    c{3}='L ';
  else
    c{3}='U ';
  end
  if trs
    a=fstdebug_sampmat(n,m);
  else
    a=fstdebug_sampmat(m,n);
  end
  u=rand;
  if u<(1/3)
    ret2=fstdeb_dsyrk(c,a,trs);
    ret1=c(:,:);
    fst_dsyrk(ret1,a,trs);
    alpha=1; beta=0;
  elseif u<(2/3)
    ret2=fstdeb_dsyrk(c,a,trs,alpha);
    ret1=c(:,:);
    fst_dsyrk(ret1,a,trs,alpha);
    beta=0;
  else
    ret2=fstdeb_dsyrk(c,a,trs,alpha,beta);
    ret1=c(:,:);
    fst_dsyrk(ret1,a,trs,alpha,beta);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dsyrk_error','c','a','trs','alpha','beta','ret1', ...
	 'ret2');
    error('ERROR (DSYRK)!');
  end
end
fprintf(1,'OK, done.\n');

% DTRMM
fprintf(1,'DTRMM tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  trs=double(rand>0.5);
  left=double(rand>0.5);
  alpha=rand;
  b=fstdebug_sampmat(m,n);
  if left
    a=fstdebug_sampmat(m,m,1);
  else
    a=fstdebug_sampmat(n,n,1);
  end
  if rand>0.5
    a{3}='L ';
  else
    a{3}='U ';
  end
  u=rand;
  if u<(1/3)
    a{3}(2)='U';
  elseif u<(2/3)
    a{3}(2)='N';
  end
  u=rand;
  if u<0.5
    ret2=fstdeb_dtrmm(b,a,left,trs);
    ret1=b(:,:);
    fst_dtrmm(ret1,a,left,trs);
    alpha=1;
  else
    ret2=fstdeb_dtrmm(b,a,left,trs,alpha);
    ret1=b(:,:);
    fst_dtrmm(ret1,a,left,trs,alpha);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dtrmm_error','b','a','trs','alpha','left','ret1', ...
	 'ret2');
    error('ERROR (DTRMM)!');
  end
end
fprintf(1,'OK, done.\n');

% DTRSM
fprintf(1,'DTRSM tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  trs=double(rand>0.5);
  left=double(rand>0.5);
  alpha=rand;
  b=fstdebug_sampmat(m,n);
  if left
    a=fstdebug_sampmat(m,m,1);
  else
    a=fstdebug_sampmat(n,n,1);
  end
  if rand>0.5
    a{3}='L ';
  else
    a{3}='U ';
  end
  u=rand;
  if u<(1/3)
    a{3}(2)='U';
  else
    % Diag. must be positive
    sz=a{2}(3);
    ys=a{2}(1); xs=a{2}(2);
    posdg=rand(sz,1)+(1e-7);
    bm=size(a{1},1);
    a{1}(((xs-1)*bm+ys):(bm+1):((xs+sz-2)*bm+ys+sz-1))=posdg;
    if u<(2/3)
      a{3}(2)='N';
    end
  end
  u=rand;
  if u<0.5
    ret2=fstdeb_dtrsm(b,a,left,trs);
    ret1=b(:,:);
    fst_dtrsm(ret1,a,left,trs);
    alpha=1;
  else
    ret2=fstdeb_dtrsm(b,a,left,trs,alpha);
    ret1=b(:,:);
    fst_dtrsm(ret1,a,left,trs,alpha);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dtrsm_error','b','a','trs','alpha','left','ret1', ...
	 'ret2');
    error('ERROR (DTRSM)!');
  end
end
fprintf(1,'OK, done.\n');

% MULDIAG
fprintf(1,'MULDIAG tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  left=double(rand>0.5);
  a=fstdebug_sampmat(m,n);
  if left
    d=2*rand(m,1)-1;
  else
    d=2*rand(n,1)-1;
  end
  ret2=fstdeb_muldiag(a,d,left);
  ret1=a(:,:);
  fst_muldiag(ret1,d,left);
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('muldiag_error','a','d','left', ...
	 'ret1','ret2');
    error('ERROR (MULDIAG)!');
  end
end
fprintf(1,'OK, done.\n');

% PERMUTE
fprintf(1,'PERMUTE tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  col=double(rand>0.5);
  a=fstdebug_sampmat(m,n);
  if col
    perm=randperm(n);
  else
    perm=randperm(m);
  end
  ret2=fstdeb_permute(a,perm,col);
  ret1=a(:,:);
  fst_permute(ret1,perm,col);
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('permute_error','a','perm','col', ...
	 'ret1','ret2');
    error('ERROR (PERMUTE)!');
  end
end
fprintf(1,'OK, done.\n');

% FLIPDIMS
fprintf(1,'FLIPDIMS tests...\n');
for i=1:num
  n=1+floor(maxn*rand);
  d1=1+floor(10*rand); d2=1+floor(10*rand);
  a=rand(n*d1,d2);
  ret2=fstdeb_flipdims(a,n);
  ret1=a(:,:);
  fst_flipdims(ret1,n);
  if max(max(reldiff(ret1,ret2,1e-15)))>(1e-7)
    save('flipdims_error','a','n','ret1','ret2');
    error('ERROR (KLR_FLIPDIMS)!');
  end
end
fprintf(1,'OK, done.\n');

% RESHAPE
fprintf(1,'RESHAPE tests...\n');
for i=1:num
  n=1+floor(maxn*rand);
  m=1+floor(maxm*rand);
  a=2*rand(m,n)-1;
  ret2=reshape(a,n,m);
  ret1=a(:,:);
  fst_reshape(ret1,n,m);
  if max(max(reldiff(ret1,ret2,1e-15)))>(1e-7)
    save('reshape_error','a','n','ret1','ret2');
    error('ERROR (KLR_RESHAPE)!');
  end
end
fprintf(1,'OK, done.\n');

% ADDVEC
fprintf(1,'ADDVEC tests...\n');
for i=1:num
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  col=double(rand>0.5);
  c=fstdebug_sampmat(m,n);
  if col
    x=2*rand(m,1)-1;
  else
    x=2*rand(n,1)-1;
  end
  alpha=rand;
  if rand<0.5
    ret2=fstdeb_addvec(c,x,col);
    ret1=c(:,:);
    fst_addvec(ret1,x,col);
    alpha=1;
  else
    ret2=fstdeb_addvec(c,x,col,alpha);
    ret1=c(:,:);
    fst_addvec(ret1,x,col,alpha);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('addvec_error','c','x','col','alpha', ...
	 'ret1','ret2');
    error('ERROR (ADDVEC)!');
  end
end
fprintf(1,'OK, done.\n');

% SUMPOS
fprintf(1,'SUMPOS tests...\n');
for i=1:num
  n=ceil(rand*maxn);
  yn=ceil(rand*15);
  x=2*rand(n,1)-1;
  ind=ceil(rand(n,1)*yn);
  ret1=fst_sumpos(x,ind,yn);
  ret2=fstdeb_sumpos(x,ind,yn);
  if max(max(reldiff(ret1,ret2,1e-15)))>(1e-7)
    save('sumpos_error','x','yn','ind', ...
	 'ret1','ret2');
    error('ERROR (SUMPOS)!');
  end
end
fprintf(1,'OK, done.\n');

% DSPAMM
fprintf(1,'DSPAMM tests...\n');
for i=1:num
  i
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  k=1+floor(rand*maxn);
  alpha=rand; beta=rand;
  tra=double(rand>0.5);
  c=fstdebug_sampmat(m,n);
  b=fstdebug_sampmat(k,n);
  if tra
    a=sprandn(k,m,1/5);
  else
    a=sprandn(m,k,1/5);
  end
  u=rand;
  if u<(1/3)
    ret2=fstdeb_dspamm(c,a,b,tra);
    ret1=c(:,:);
    fst_dspamm(ret1,a,b,tra);
    alpha=1; beta=0;
  elseif u<(2/3)
    ret2=fstdeb_dspamm(c,a,b,tra,alpha);
    ret1=c(:,:);
    fst_dspamm(ret1,a,b,tra,alpha);
    beta=0;
  else
    ret2=fstdeb_dspamm(c,a,b,tra,alpha,beta);
    ret1=c(:,:);
    fst_dspamm(ret1,a,b,tra,alpha,beta);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dspamm_error','c','a','b','tra','alpha','beta', ...
	 'ret1','ret2');
    error('ERROR (DSPAMM)!');
  end
end

% DSPAMMT
fprintf(1,'DSPAMMT tests...\n');
for i=1:num
  i
  m=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  k=1+floor(rand*maxn);
  alpha=rand; beta=rand;
  c=fstdebug_sampmat(m,k);
  b=fstdebug_sampmat(m,k);
  a=sprandn(m,n,1/5);
  u=rand;
  if u<(1/3)
    ret2=fstdeb_dspammt(c,a,b);
    ret1=c(:,:);
    fst_dspammt(ret1,a,b);
    alpha=1; beta=0;
  elseif u<(2/3)
    ret2=fstdeb_dspammt(c,a,b,alpha);
    ret1=c(:,:);
    fst_dspammt(ret1,a,b,alpha);
    beta=0;
  else
    ret2=fstdeb_dspammt(c,a,b,alpha,beta);
    ret1=c(:,:);
    fst_dspammt(ret1,a,b,alpha,beta);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dspammt_error','c','a','b','tra','alpha','beta', ...
	 'ret1','ret2');
    error('ERROR (DSPAMMT)!');
  end
end

% DSPAMMT2
fprintf(1,'DSPAMMT2 tests...\n');
for i=1:num
  i
  m1=1+floor(rand*maxm);
  m2=1+floor(rand*maxm);
  n=1+floor(rand*maxn);
  k=1+floor(rand*maxn);
  alpha=rand; beta=rand;
  d=fstdebug_sampmat(m1,k);
  c=fstdebug_sampmat(m2,k);
  a=sprandn(m1,n,1/5);
  b=sprandn(m2,n,1/8);
  u=rand;
  if u<(1/3)
    ret2=fstdeb_dspammt2(d,a,b,c);
    ret1=d(:,:);
    fst_dspammt2(ret1,a,b,c);
    alpha=1; beta=0;
  elseif u<(2/3)
    ret2=fstdeb_dspammt2(d,a,b,c,alpha);
    ret1=d(:,:);
    fst_dspammt2(ret1,a,b,c,alpha);
    beta=0;
  else
    ret2=fstdeb_dspammt2(d,a,b,c,alpha,beta);
    ret1=d(:,:);
    fst_dspammt2(ret1,a,b,c,alpha,beta);
  end
  if max(max(reldiff(fstdebug_getbuffmat(ret1),ret2,1e-15)))>(1e-7)
    save('dspammt2_error','d','c','a','b','tra','alpha','beta', ...
	 'ret1','ret2');
    error('ERROR (DSPAMMT2)!');
  end
end

% ACCUMULATE
fprintf(1,'ACCUMULATE tests...\n');
for i=1:num
  i
  n=1+floor(rand*maxn);
  nc=1+floor(rand*30);
  a=fstdebug_sampmat(n,nc);
  nb=1+floor(rand*30);
  b=[];
  for j=1:nb
    sz=ceil(rand*10);
    b{j}=ceil(rand(sz,1)*nc);
  end
  iv=ceil(rand(n,1)*nb);
  ret1=fst_accumulate(a,iv,b);
  ret2=fstdeb_accumulate(a,iv,b);
  if max(max(reldiff(ret1,ret2,1e-15)))>(1e-7)
    save('accumulate_error','a','iv','b','ret1','ret2');
    error('ERROR (ACCUMULATE)!');
  end
end
