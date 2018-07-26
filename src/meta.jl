


function create_quote(fastmath = false, inbounds = true)
    if fastmath == true && inbounds == true
        q = quote
            @fastmath @inbounds begin
            end
        end
        qa = q.args[2].args[3].args[3].args
    elseif inbounds == true
        q = quote
            @inbounds begin
            end
        end
        qa = q.args[2].args[3].args
    elseif fastmath == true
        q = quote
            @fastmath begin
            end
        end
        qa = q.args[2].args[3].args
    else
        q = quote end
        qa = q.args
    end
    q, qa
end

inv_expr(x::Union{Symbol,Expr}) = :(inv($x))
inv_expr(x) = inv(x)
exponentiate_expr(x::Union{Symbol,Expr}, y) = :($x^$y)
exponentiate_expr(x, y::Union{Symbol,Expr}) = :($x^$y)
exponentiate_expr(x::Union{Symbol,Expr}, y::Union{Symbol,Expr}) = :($x^$y)
exponentiate_expr(x, y) = x^y

"""
Names are mangled by appending ndq, in case you want to access them within the body.
Arguments:
Body is the body expression
K is the degree
term is the symbol for the polynomial term.
D is the maximum degree of the polynomial.
L is the Lebesque measure of the polynomial.

Note that the terms will be called
term_k, k = 1,...,K
so you can make use of that in the body.

Note: these generated functions are a hack.
They really ought to be replaced by implementing the iteration interface.
"""
function norm_degree_quote(body, K, term = :term, D = :D, L = :L)
    Kp1 = K+1
    quote
        Linvndq = $(inv_expr(L))
        jndq_0 = $D
        DLndq = $(exponentiate_expr(D, L))
        sumndq_0 = 0
        @nloops $K $term kndq -> begin
            0:jndq_{$K-kndq}
        end kndq -> begin
            sumndq_{$Kp1-kndq} = sumndq_{$K-kndq} + $(Symbol(term, :_kndq))^$L
            jndq_{$Kp1-kndq} = floor(Int, (DLndq - sumndq_{$Kp1-kndq})^Linvndq)
        end begin
            $body
        end
    end
end

function max_degree_quote(body, K, term = :term, D = :D)
    quote
        rmdq_0 = $D
        Base.Cartesian.@nloops $K $term kmdq-> begin
            rmdq_{$K-kmdq}:-1:0
        end kmdq -> begin
            rmdq_{$(K+1)-j} = rmdq_{$N-j} - $(Symbol(term, :_kmdq))
        end begin
            $body
        end
    end
end

"""
Takes a function that gets passed `sym` and `n` as two argument to choose the iteration range.
"""
function fbound_nloops(N, sym, f)
    q, qa = create_quote(false, false)
    for n ∈ 1:N
        push!(qa, :(for $(Symbol(sym, :_, n)) in $(f(sym, n)); end ))
        qa = qa[2].args[2].args
    end
    q, qa
end

function constrained_nloops(N, sym, constraint, lb = 0)
    q, qa = create_quote(false, false)
    push!(qa, :($(Symbol(sym, :_itersumremainder_, 0)) = $constraint))
    for n ∈ 1:N
        push!(qa, :(for $(Symbol(sym, :_, n)) in $(lb):$(Symbol(sym, :_itersumremainder_, n-1)); end ))
        # @show length(qa), qa[3]
        # @show qa[3].args[2]
        qa = qa[3].args[2].args
        push!(qa, :($(Symbol(sym, :_itersumremainder_, n)) = $(Symbol(sym, :_itersumremainder_, n-1)) - $(Symbol(sym, :_, n)) ))
    end
    q, qa
end
