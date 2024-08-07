double Legendre_polynomials(int l, double x)
{
    switch(l)
    {
        case 0 :
            return 1;
        break;

        case 1 : 
            return x;
        break;

        case 2 :
            return 0.5*(3.0*x*x-1.0);
        break;

        case 3 :
            return 0.5*(5.0*x*x*x-3.0*x);
        break;

        case 4 :
            return 0.125*(35.0*x*x*x*x-30.0*x*x+3.0);
        break;
    }
    
    return 0;
}