function f1(x)
    return sum([1,-1].*x.*x)
end

function f2(x)
    return sum(x.*x)
end

function f3(x)
    return sum([1,0].*x)*sum([0,1].*x)
end

function f4(x)
    return sum([1,0].*x)*exp(-sum(abs2,x))
end

function f5(x)
    return sum([1,0].*x)*sin(sum([0,1].*x))
end

function f6(x)
    return sin(sum([1,0].*x))*sum(abs2,[0,1].*x)/sum([1,0].*x)
end

function f7(x)
    return sum([1,0].*x)*sin(sum([1,0].*x))*sin(sum([0,1].*x))
end

function f8(x)
    return sum([1,-1].*x.*x)^2
end

function f9(x)
    return exp(-10*sum([1,-1].*x.*x)^2)
end

function f0(x)
    return sqrt(sum(x.*x))
end

function f_inverse(x)
    return sum([1,0].*x)/sum([0,1].*x)
end
