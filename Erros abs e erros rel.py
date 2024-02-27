import math
x = 1
x_ant = x
x_exa = (1+math.sqrt(5))/2
    
for k in range(20):
        x = 1+1/x
        dif = abs(x-x_ant)/abs(x)
        err = abs(x-x_exa)/abs(x_exa)
        print (x, x_exa, err, dif)
        x_ant = x
