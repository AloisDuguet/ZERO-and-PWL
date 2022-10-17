function factorial(x)
    if x==0
        return 1
    else
        return x * factorial(x-=1)
    end
end

function DLcos(x,ordre)
    if ordre == 0
        return 1
    else
        return (-1)^(ordre) * x^(2*ordre)/factorial(2*ordre) + DLcos(x,ordre-=1)
    end
end

function DLsin(x,ordre)
    if ordre == 0
        return x
    else
        return (-1)^(ordre) * x^(2*ordre+1)/factorial(2*ordre+1) + DLsin(x,ordre-=1)
    end
end

function ftest1(x,y)
	x^2 - y^2
end

function ftest2(x,y)
	return x^2 + y^2
end

function ftest3(x,y)
	return x * y
end

function ftest4(x,y)
	return x * exp(-x^2-y^2)
end

function ftest5(x,y)
	return x * (y - y^3/6 + y^5/120 - y^7/5040 + y^9/362880 - y^11/39916800 + y^13/6227020800 - y^15/1307674368000 + y^17/355687428096000 - y^19/121645100408832000)
end

function ftest6(x,y)
	return (x - x^3/6 + x^5/120 - x^7/5040 + x^9/362880 - x^11/39916800 + x^13/6227020800 - x^15/1307674368000 + x^17/355687428096000 - x^19/121645100408832000)/x * y^2
end

function ftest7(x,y)
	return x * (x - x^3/6 + x^5/120 - x^7/5040 + x^9/362880 - x^11/39916800 + x^13/6227020800 - x^15/1307674368000 + x^17/355687428096000 - x^19/121645100408832000) * (y - y^3/6 + y^5/120 - y^7/5040 + y^9/362880 - y^11/39916800 + y^13/6227020800 - y^15/1307674368000 + y^17/355687428096000 - y^19/121645100408832000)
end

function ftest8(x,y)
	return (x^2 - y^2)^2
end

function ftest9(x,y)
	return exp(-10 * (x^2 - y^2)^2)
end

function ftest10(x,y)
    return sqrt(sum(x^2+y^2))
end
