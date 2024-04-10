from Crypto.Util.number import *
def coron(pol, X, Y, k=2, debug=False):
    """
    Returns all small roots of pol.

    Applies Coron's reformulation of Coppersmith's algorithm for finding small
    integer roots of bivariate polynomials modulo an integer.

    Args:
        pol: The polynomial to find small integer roots of.
        X: Upper limit on x.
        Y: Upper limit on y.
        k: Determines size of lattice. Increase if the algorithm fails.
        debug: Turn on for debug print stuff.

    Returns:
        A list of successfully found roots [(x0,y0), ...].

    Raises:
        ValueError: If pol is not bivariate
    """

    if pol.nvariables() != 2:
        raise ValueError("pol is not bivariate")

    P.<x,y> = PolynomialRing(ZZ)
    pol = pol(x,y)

    # Handle case where pol(0,0) == 0
    xoffset = 0

    while pol(xoffset,0) == 0:
        xoffset += 1

    pol = pol(x+xoffset,y)

    # Handle case where gcd(pol(0,0),X*Y) != 1
    while gcd(pol(0,0), X) != 1:
        X = next_prime(X, proof=False)

    while gcd(pol(0,0), Y) != 1:
        Y = next_prime(Y, proof=False)

    pol = P(pol/gcd(pol.coefficients())) # seems to be helpful
    p00 = pol(0,0)
    delta = max(pol.degree(x),pol.degree(y)) # maximum degree of any variable

    W = max(abs(i) for i in pol(x*X,y*Y).coefficients())
    u = W + ((1-W) % abs(p00))
    N = u*(X*Y)^k # modulus for polynomials

    # Construct polynomials
    p00inv = inverse_mod(p00,N)
    polq = P(sum((i*p00inv % N)*j for i,j in zip(pol.coefficients(),
                                                 pol.monomials())))
    polynomials = []
    for i in range(delta+k+1):
        for j in range(delta+k+1):
            if 0 <= i <= k and 0 <= j <= k:
                polynomials.append(polq * x^i * y^j * X^(k-i) * Y^(k-j))
            else:
                polynomials.append(x^i * y^j * N)

    # Make list of monomials for matrix indices
    monomials = []
    for i in polynomials:
        for j in i.monomials():
            if j not in monomials:
                monomials.append(j)
    monomials.sort()

    # Construct lattice spanned by polynomials with xX and yY
    L = matrix(ZZ,len(monomials))
    for i in range(len(monomials)):
        for j in range(len(monomials)):
            L[i,j] = polynomials[i](X*x,Y*y).monomial_coefficient(monomials[j])

    # makes lattice upper triangular
    # probably not needed, but it makes debug output pretty
    L = matrix(ZZ,sorted(L,reverse=True))

    if debug:
        print("Bitlengths of matrix elements (before reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

    L = L.LLL()

    if debug:
        print("Bitlengths of matrix elements (after reduction):")
        print(L.apply_map(lambda x: x.nbits()).str())

    roots = []

    for i in range(L.nrows()):
        if debug:
            print("Trying row {}".format(i))

        # i'th row converted to polynomial dividing out X and Y
        pol2 = P(sum(map(mul, zip(L[i],monomials)))(x/X,y/Y))

        r = pol.resultant(pol2, y)

        if r.is_constant(): # not independent
            continue

        for x0, _ in r.univariate_polynomial().roots():
            if x0-xoffset in [i[0] for i in roots]:
                continue
            if debug:
                print("Potential x0:",x0)
            for y0, _ in pol(x0,y).univariate_polynomial().roots():
                if debug:
                    print("Potential y0:",y0)
                if (x0-xoffset,y0) not in roots and pol(x0,y0) == 0:
                    roots.append((x0-xoffset,y0))
    return roots

def main():
    N  = 138963980427736364646203557164328211078134463518489686365728312873583832517087170768576679472472907142081360480944201759920246566585465801088226164314480607014663211599932950864391702460227584467326051919881067028851940610382044445003060103566003934601979805899293539507221062915314813557293919231917284247667
    n  = 1485715964481761497309522733620825737885569961284688766942216863704985393094065876545992131370884059645617234469978112000000000000000000000
    c = 26363325527372681448374836719361674028908733933823971039273016094221739663363697355984980560218941405351917768372297139270315950803631724328547161889191685480725185971092638691575587334307068143724069148715129866085595445974433311000459043513392513632639058879350662222598941781017396217632160254074487773693
    a = N % n
    X = Y = 2^53 # bounds on x0 and y0

    P.<x,y> = PolynomialRing(ZZ)
    pol = (n * x + 1)*(n * y + a) - N # Should have a root at (x0,y0)

    x0, y0 = coron(pol, X, Y, k=2, debug=True)[0]
    p = n * x0 + 1
    q = n * y0 + a
    print('Recovered:')
    print('x0 =',x0)
    print('y0 =',y0)
    d = pow(65537, -1, N - int(p) - int(q) + 1)
    flag = pow(c, d, N)
    print(long_to_bytes(int(flag)).decode())


if __name__ == '__main__':
    main()
