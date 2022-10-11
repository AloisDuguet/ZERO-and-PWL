using LinA

function find_1d_domain_to_refine(pwl, sol_val)
    # return t1, t2 and pieces_number.
    # describing the domain [t1,t2] of the pieces containing sol_val and the numbers of the corresponding pieces in list pieces_number

    return t1,t2,pieces_number
end

function refine_pwl1d(pwl, exprf, sol_val, err)
    # return pwl of one variable replaced by a piecewise linear function satisfying error err for function with expression exprf on the domain containing all pieces containing sol_val

    return pwl
end

function refine_pwl2d(pwl, str_exprf, sol_val, err)
    # return pwl of two variables replaced by a piecewise linear function satisfying error err for function with expression exprf on the domain containing all pieces containing sol_val

    return pwl
end
