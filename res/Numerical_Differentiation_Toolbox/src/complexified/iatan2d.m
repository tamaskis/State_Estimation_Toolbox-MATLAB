% https://mdolab.engin.umich.edu/misc/files/complexify.f90
function result = iatan2d(y,x)
    a = real(y);
    b = imag(y);
    c = real(x);
    d = imag(x);
    result = (180/pi)*atan2(a,c)+1i*((c*b-a*d)/(a^2+c^2));
end