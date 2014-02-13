cp lapacke.h a
indent -kr -l50000 a
grep LAPACKE_d a > b
mv b a
grep -v _work a > b
mv b a
grep -v float a > b
mv b a
grep -v lapack_logical a > b
mv b a
sed -i 's/lapack_int/int/g' a
indent -kr -l50000 a
sed -i '/^$/d' a
grep -v SELECT a > b
mv b a
sed -i 's/LAPACKE_/mql_/g' a
sed -i 's/const //g' a
sed -i 's/double \*a,/double \&a[],/g' a
sed -i 's/double \*a)/double \&a[])/g' a
sed -i 's/double \*b,/double \&b[],/g' a
sed -i 's/double \*b)/double \&b[])/g' a

sed -i 's/double \*bb,/double \&bb[],/g' a
sed -i 's/double \*bb)/double \&bb[])/g' a

sed -i 's/double \*c,/double \&c[],/g' a
sed -i 's/double \*c)/double \&c[])/g' a

sed -i 's/double \*d,/double \&d[],/g' a
sed -i 's/double \*d)/double \&d[])/g' a

sed -i 's/double \*e,/double \&e[],/g' a
sed -i 's/double \*e)/double \&e[])/g' a

sed -i 's/double \*f,/double \&f[],/g' a
sed -i 's/double \*f)/double \&f[])/g' a

sed -i 's/double \*p,/double \&p[],/g' a
sed -i 's/double \*p)/double \&p[])/g' a

sed -i 's/double \*q,/double \&q[],/g' a
sed -i 's/double \*q)/double \&q[])/g' a

sed -i 's/double \*r,/double \&r[],/g' a
sed -i 's/double \*r)/double \&r[])/g' a

sed -i 's/double \*s,/double \&s[],/g' a
sed -i 's/double \*s)/double \&s[])/g' a

sed -i 's/double \*h,/double \&h[],/g' a
sed -i 's/double \*h)/double \&h[])/g' a

sed -i 's/double \*ef,/double \&ef[],/g' a
sed -i 's/double \*ef)/double \&ef[])/g' a

sed -i 's/double \*t,/double \&t[],/g' a
sed -i 's/double \*t)/double \&t[])/g' a
sed -i 's/double \*u,/double \&u[],/g' a
sed -i 's/double \*u)/double \&u[])/g' a
sed -i 's/double \*v,/double \&v[],/g' a
sed -i 's/double \*v)/double \&v[])/g' a
sed -i 's/double \*w,/double \&w[],/g' a
sed -i 's/double \*w)/double \&w[])/g' a
sed -i 's/double \*x,/double \&x[],/g' a
sed -i 's/double \*x)/double \&x[])/g' a
sed -i 's/double \*y,/double \&y[],/g' a
sed -i 's/double \*y)/double \&y[])/g' a
sed -i 's/double \*z,/double \&z[],/g' a
sed -i 's/double \*z)/double \&z[])/g' a
sed -i 's/double \*vt,/double \&vt[],/g' a
sed -i 's/double \*vt)/double \&vt[])/g' a
sed -i 's/double \*tau,/double \&tau[],/g' a
sed -i 's/double \*tau)/double \&tau[])/g' a

sed -i 's/double \*taua,/double \&taua[],/g' a
sed -i 's/double \*taua)/double \&taua[])/g' a

sed -i 's/double \*taub,/double \&taub[],/g' a
sed -i 's/double \*taub)/double \&taub[])/g' a

sed -i 's/double \*taup,/double \&taup[],/g' a
sed -i 's/double \*taup)/double \&taup[])/g' a

sed -i 's/double \*tauq,/double \&tauq[],/g' a
sed -i 's/double \*tauq)/double \&tauq[])/g' a

sed -i 's/double \*ab,/double \&ab[],/g' a
sed -i 's/double \*ab)/double \&ab[])/g' a
sed -i 's/double \*af,/double \&af[],/g' a
sed -i 's/double \*af)/double \&af[])/g' a
sed -i 's/int \*ipiv)/int \&ipiv[])/g' a
sed -i 's/int \*piv/int \&piv[]/g' a
sed -i 's/int \*iq)/int \&iq[])/g' a
sed -i 's/int \*ipiv,/int \&ipiv[],/g' a
sed -i 's/double \*ap,/double \&ap[],/g' a
sed -i 's/double \*ap)/double \&ap[])/g' a

sed -i 's/double \*afp,/double \&afp[],/g' a
sed -i 's/double \*afp)/double \&afp[])/g' a

sed -i 's/double \*bp,/double \&bp[],/g' a
sed -i 's/double \*bp)/double \&bp[])/g' a

sed -i 's/double \*afb,/double \&afb[],/g' a
sed -i 's/double \*afb)/double \&afb[])/g' a

sed -i 's/int \*jpvt,/int \&jpvt[],/g' a
sed -i 's/int \*jpvt)/int \&jpvt[])/g' a

sed -i 's/int \*ifail,/int \&ifail[],/g' a
sed -i 's/int \*ifail)/int \&ifail[])/g' a

sed -i 's/double \*arf,/double \&arf[],/g' a
sed -i 's/double \*arf)/double \&arf[])/g' a

sed -i 's/double \*sep,/double \&sep[],/g' a
sed -i 's/double \*sep)/double \&sep[])/g' a

sed -i 's/double \*rcond/double \&rcond/g' a
sed -i 's/double \*scond/double \&scond/g' a

sed -i 's/double \*rowcnd/double \&rowcnd/g' a
sed -i 's/double \*colcnd/double \&colcnd/g' a

sed -i 's/double \*lscale/double \&lscale/g' a
sed -i 's/double \*rscale/double \&rscale/g' a
sed -i 's/double \*rpvgrw/double \&rpvgrw/g' a

sed -i 's/double \*wi,/double \&wi[],/g' a
sed -i 's/double \*wi)/double \&wi[])/g' a

sed -i 's/double \*wr,/double \&wr[],/g' a
sed -i 's/double \*wr)/double \&wr[])/g' a

sed -i 's/double \*vr,/double \&vr[],/g' a
sed -i 's/double \*vr)/double \&vr[])/g' a
sed -i 's/double \*vl,/double \&vl[],/g' a
sed -i 's/double \*vl)/double \&vl[])/g' a

sed -i 's/double \*pt,/double \&pt[],/g' a
sed -i 's/double \*pt)/double \&pt[])/g' a

sed -i 's/double \*scale,/double \&scale[],/g' a
sed -i 's/double \*scale)/double \&scale[])/g' a

sed -i 's/double \*params,/double \&params[],/g' a
sed -i 's/double \*params)/double \&params[])/g' a

sed -i 's/double \*alpha,/double \&alpha[],/g' a
sed -i 's/double \*alpha)/double \&alpha[])/g' a

sed -i 's/double \*alphai,/double \&alphai[],/g' a
sed -i 's/double \*alphai)/double \&alphai[])/g' a

sed -i 's/double \*alphar,/double \&alphar[],/g' a
sed -i 's/double \*alphar)/double \&alphar[])/g' a

sed -i 's/double \*beta,/double \&beta[],/g' a
sed -i 's/double \*beta)/double \&beta[])/g' a

sed -i 's/double \*x12,/double \&x12[],/g' a
sed -i 's/double \*x12)/double \&x12[])/g' a

sed -i 's/double \*x11,/double \&x11[],/g' a
sed -i 's/double \*x11)/double \&x11[])/g' a

sed -i 's/double \*x21,/double \&x21[],/g' a
sed -i 's/double \*x21)/double \&x21[])/g' a

sed -i 's/double \*x22,/double \&x22[],/g' a
sed -i 's/double \*x22)/double \&x22[])/g' a

sed -i 's/double \*u2,/double \&u2[],/g' a
sed -i 's/double \*u2)/double \&u2[])/g' a

sed -i 's/double \*u1,/double \&u1[],/g' a
sed -i 's/double \*u1)/double \&u1[])/g' a

sed -i 's/double \*v1t,/double \&v1t[],/g' a
sed -i 's/double \*v1t)/double \&v1t[])/g' a

sed -i 's/double \*v2t,/double \&v2t[],/g' a
sed -i 's/double \*v2t)/double \&v2t[])/g' a

sed -i 's/int \*iseed,/int \&iseed[],/g' a
sed -i 's/int \*iseed)/int \&iseed[])/g' a

sed -i 's/int \*ncycle,/int \&ncycle[],/g' a
sed -i 's/int \*ncycle)/int \&ncycle[])/g' a

sed -i 's/int \*iter,/int \&iter[],/g' a
sed -i 's/int \*iter)/int \&iter[])/g' a

sed -i 's/double \*amax/double \&amax/g' a

sed -i 's/int \*rank/int \&rank/g' a

sed -i 's/int \*ilo/int \&ilo/g' a
sed -i 's/int \*ihi/int \&ihi/g' a

sed -i 's/int \*ifst/int \&ifst/g' a
sed -i 's/int \*ilst/int \&ilst/g' a

sed -i 's/int \*isuppz/int \&isuppz/g' a

sed -i 's/int \*iblock/int \&iblock[]/g' a
sed -i 's/int \*isplit/int \&isplit[]/g' a

sed -i 's/int \*nsplit/int \&nsplit/g' a
sed -i 's/int \*ifailv/int \&ifailv[]/g' a

sed -i 's/int \*m,/int \&m,/g' a
sed -i 's/int \*k,/int \&k,/g' a
sed -i 's/int \*l,/int \&l,/g' a

sed -i 's/char \*/char /g' a

sed -i 's/double \*ferr,/double \&ferr[],/g' a
sed -i 's/double \*ferr)/double \&ferr[])/g' a

sed -i 's/double \*berr,/double \&berr[],/g' a
sed -i 's/double \*berr)/double \&berr[])/g' a

sed -i 's/double \*superb,/double \&superb[],/g' a
sed -i 's/double \*superb)/double \&superb[])/g' a

sed -i 's/double \*err_bnds_norm,/double \&err_bnds_norm[],/g' a
sed -i 's/double \*err_bnds_norm)/double \&err_bnds_norm[])/g' a

sed -i 's/double \*err_bnds_comp,/double \&err_bnds_comp[],/g' a
sed -i 's/double \*err_bnds_comp)/double \&err_bnds_comp[])/g' a

sed -i 's/double \*sva,/double \&sva[],/g' a
sed -i 's/double \*sva)/double \&sva[])/g' a

sed -i 's/double \*cs,/double \&cs,/g' a
sed -i 's/double \*cs)/double \&cs)/g' a

sed -i 's/double \*sn,/double \&sn,/g' a
sed -i 's/double \*sn)/double \&sn)/g' a

sed -i 's/double \*taup1,/double \&taup1[],/g' a
sed -i 's/double \*taup1)/double \&taup1[])/g' a
sed -i 's/double \*taup2,/double \&taup2[],/g' a
sed -i 's/double \*taup2)/double \&taup2[])/g' a

sed -i 's/double \*tauq1,/double \&tauq1[],/g' a
sed -i 's/double \*tauq1)/double \&tauq1[])/g' a
sed -i 's/double \*tauq2,/double \&tauq2[],/g' a
sed -i 's/double \*tauq2)/double \&tauq2[])/g' a

sed -i 's/double \*theta,/double \&theta[],/g' a
sed -i 's/double \*theta)/double \&theta[])/g' a

sed -i 's/double \*phi,/double \&phi[],/g' a
sed -i 's/double \*phi)/double \&phi[])/g' a

sed -i 's/double \*df,/double \&df[],/g' a
sed -i 's/double \*df)/double \&df[])/g' a

sed -i 's/double \*dl,/double \&dl[],/g' a
sed -i 's/double \*dl)/double \&dl[])/g' a

sed -i 's/double \*dlf,/double \&dlf[],/g' a
sed -i 's/double \*dlf)/double \&dlf[])/g' a

sed -i 's/double \*du,/double \&du[],/g' a
sed -i 's/double \*du)/double \&du[])/g' a

sed -i 's/double \*du2,/double \&du2[],/g' a
sed -i 's/double \*du2)/double \&du2[])/g' a

sed -i 's/double \*duf,/double \&duf[],/g' a
sed -i 's/double \*duf)/double \&duf[])/g' a

sed -i 's/double \*b11d,/double \&b11d[],/g' a
sed -i 's/double \*b11d)/double \&b11d[])/g' a

sed -i 's/double \*b12d,/double \&b12d[],/g' a
sed -i 's/double \*b12d)/double \&b12d[])/g' a

sed -i 's/double \*b22d,/double \&b22d[],/g' a
sed -i 's/double \*b22d)/double \&b22d[])/g' a

sed -i 's/double \*b22e,/double \&b22e[],/g' a
sed -i 's/double \*b22e)/double \&b22e[])/g' a

sed -i 's/double \*b11e,/double \&b11e[],/g' a
sed -i 's/double \*b11e)/double \&b11e[])/g' a

sed -i 's/double \*b12e,/double \&b12e[],/g' a
sed -i 's/double \*b12e)/double \&b12e[])/g' a

sed -i 's/double \*b21d,/double \&b21d[],/g' a
sed -i 's/double \*b21d)/double \&b21d[])/g' a

sed -i 's/double \*b21e,/double \&b21e[],/g' a
sed -i 's/double \*b21e)/double \&b21e[])/g' a

sed -i 's/double \*dif/double \&dif/g' a

sed -i 's/double \*rpivot,/double \&rpivot[],/g' a
sed -i 's/double \*rpivot)/double \&rpivot[])/g' a

sed -i 's/double \*stat,/double \&stat[],/g' a
sed -i 's/double \*stat)/double \&stat[])/g' a

sed -i 's/int \*istat,/int \&istat[],/g' a
sed -i 's/int \*istat)/int \&istat[])/g' a

sed -i 's/int \*isgn,/int \&isgn[],/g' a
sed -i 's/int \*isgn)/int \&isgn[])/g' a

sed -i 's/int \*isave,/int \&isave[],/g' a
sed -i 's/int \*isave)/int \&isave[])/g' a

sed -i 's/int \*kase,/int \&kase,/g' a
sed -i 's/int \*kase)/int \&kase)/g' a

sed -i 's/int \*iwork,/int \&iwork[],/g' a
sed -i 's/int \*iwork)/int \&iwork[])/g' a

sed -i 's/double \*work,/double \&work[],/g' a
sed -i 's/double \*work)/double \&work[])/g' a

sed -i 's/double \*abnrm/double \&abnrm/g' a
sed -i 's/double \*bbnrm/double \&bbnrm/g' a

sed -i 's/double \*est/double \&est/g' a

cp a aa
sed -i 's/int matrix_order, //g' a

indent -kr -l50000 a
indent -kr -l50000 aa
sed -i '/^$/d' a
sed -i '/^$/d' aa

echo "non    :"; grep \* a | wc -l
echo "success:"; grep -v \* a | wc -l

sed -i 's/int mql_/\treturn LAPACKE_/g' aa
sed -i 's/double mql_/\treturn LAPACKE_/g' aa
sed -i 's/char //g' aa
sed -i 's/int n,/(lapack_int) n,/g' aa
sed -i 's/int m,/(lapack_int) m,/g' aa
sed -i 's/int k,/(lapack_int) k,/g' aa
sed -i 's/int l,/(lapack_int) l,/g' aa
sed -i 's/int lda,/(lapack_int) lda,/g' aa
sed -i 's/int ldab,/(lapack_int) ldab,/g' aa
sed -i 's/int ldab)/(lapack_int) ldab)/g' aa
sed -i 's/int ldb,/(lapack_int) ldb,/g' aa
sed -i 's/int ldb)/(lapack_int) ldb)/g' aa
sed -i 's/int ldc,/(lapack_int) ldc,/g' aa
sed -i 's/int ldc)/(lapack_int) ldc)/g' aa
sed -i 's/int ku,/(lapack_int) ku,/g' aa
sed -i 's/int kl,/(lapack_int) kl,/g' aa
sed -i 's/int ldu,/(lapack_int) ldu,/g' aa
sed -i 's/int ldv,/(lapack_int) ldv,/g' aa
sed -i 's/int ldvt,/(lapack_int) ldvt,/g' aa
sed -i 's/int ldt,/(lapack_int) ldt,/g' aa
sed -i 's/int ldt)/(lapack_int) ldt)/g' aa
sed -i 's/int ldq,/(lapack_int) ldq,/g' aa
sed -i 's/int matrix_order,/LAPACK_COL_MAJOR,/g' aa
sed -i 's/double \&a\[\]/a/g' aa
sed -i 's/double \&ab\[\]/ab/g' aa
sed -i 's/double \&b\[\]/b/g' aa
sed -i 's/double \&c\[\]/c/g' aa
sed -i 's/double \&d\[\]/d/g' aa
sed -i 's/double \&e\[\]/e/g' aa
sed -i 's/double \&u\[\]/u/g' aa
sed -i 's/double \&q\[\]/q/g' aa
sed -i 's/double \&v\[\]/v/g' aa
sed -i 's/double \&vt\[\]/vt/g' aa
sed -i 's/double \&t\[\]/t/g' aa
sed -i 's/double \&x\[\]/x/g' aa
sed -i 's/double \&ap\[\]/ap/g' aa
sed -i 's/double \&//g' aa
sed -i 's/int \&//g' aa
sed -i 's/int /(lapack_int) /g' aa
sed -i 's/\[\]//g' aa

sed -i 's/;/;\n/g'  aa

sed -i 's/^int/_DLLAPI int __stdcall/g' a 
sed -i 's/^double/_DLLAPI double __stdcall/g' a
sed -i 's/;/\n{\n}\n/g'  a
