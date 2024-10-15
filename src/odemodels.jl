
function irmodel!(du, u, p ,t)
    I, R = u 
    γ = p[1]

    du[1] = -γ * I 
    du[2] = γ * I 
end 

function iirmodel!(du, u, p ,t)
    I1, I2, R = u 
    γ = p[1]

    du[1] = -2γ * I1 
    du[2] = 2γ * (I1 - I2) 
    du[3] = 2γ * I2 
end 
