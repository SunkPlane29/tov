# integrates the given function from a to b using the trapezoid method
# this possibly introduces some considerable error
# TODO: make a function with a more sophisticated method
function integrate_trap(f::Function, a::Real, b::Real, n::Real=100000000)::Real
    Δx = (b-a)/n
    integral = f(a) + f(b)

    for i = 1:n-1
        integral += 2f(a + i*Δx)
    end

    integral *= Δx/2

    return integral
end
