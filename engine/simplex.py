from copy import deepcopy
from numpy import argmin
from matrix import transpose, gauss, multiply, column

class ProblemType:
    MIN = -1
    MAX = 1
    
    enum = {'min' : MIN, 'max' : MAX}
    enum_inv = dict(zip(enum.values(), enum.keys()))


class Operator:
    EQUAL = 0
    LESS = -1
    GREATER = 1
    LESS_EQUAL = -2
    GREATER_EQUAL = 2
    
    enum = {'=' : EQUAL, '<' : LESS, '>' : GREATER, '<=' : LESS_EQUAL, '>=' : GREATER_EQUAL}
    enum_inv = dict(zip(enum.values(), enum.keys()))


class LinearProblem:
    
    def __init__(self, ptype, c, a, b, operators):
        self.ptype = ptype
        self.c = c
        self.a = a
        self.b = b
        self.num_real_vars = len(c)
        self.num_slack_vars = 0
        self.num_artificial_vars = 0
        self.operators = operators
        self.is_artificial = False
        self.phase_one = False
        self.is_max = False  
        self._standard_form()
    
    def _standard_form(self):        
        if self.ptype is ProblemType.MAX:
            self.ptype, self.c = -self.ptype, map(lambda(x): -x, self.c)
            self.is_max = True

        count_equals = self.operators.count(Operator.EQUAL)
        j = 0
        
        for i, op in enumerate(self.operators):
            if self.b[i] < 0.0:
                self.b[i], self.a[i], op = -self.b[i], map(lambda(x): -x, self.a[i]), -op
            elif self.operators[i] >= Operator.EQUAL:
                self.phase_one = True
            
            slack_vars = [0.0] * (len(self.operators) - count_equals)
            
            if op is not Operator.EQUAL:
                if op is Operator.LESS or op is Operator.LESS_EQUAL:
                    slack_vars[j] = 1.0
                else:
                    slack_vars[j] = -1.0
                j += 1

            self.a[i] += slack_vars
            self.operators[i] = Operator.EQUAL
                        
        self.c += [0.0] * (len(self.a) - count_equals)
        self.num_slack_vars = (len(self.operators) - count_equals) 
    
    def __str__(self):
        c = str(self.c)
        rest = '\n'.join(str(self.a[i]) + ' ' +\
                         Operator.enum_inv[self.operators[i]] + ' ' +\
                         str(self.b[i]) for i in range(len(self.a)))

        return ProblemType.enum_inv[self.ptype] + ' ' + c + ' subject to:\n' + rest
    
    @classmethod
    def artificial_problem(cls, lp):
        if lp.phase_one and not lp.is_artificial:
            artificial_lp = deepcopy(lp)
            artificial_lp.c = [0.0] * len(lp.c) + [1.0] * len(lp.a)
            artificial_lp.num_artificial_vars = len(artificial_lp.a)
            
            for i in range(len(artificial_lp.a)):
                artificial_vars = [0.0] * len(artificial_lp.a)
                artificial_vars[i] = 1.0
                artificial_lp.a[i] += artificial_vars

            artificial_lp.is_artificial = True
            return artificial_lp
        return None


def simplex(lp):
    if lp.phase_one:
        artificial_lp = LinearProblem.artificial_problem(lp)
        begin = artificial_lp.num_real_vars + artificial_lp.num_slack_vars
        end = artificial_lp.num_real_vars + artificial_lp.num_slack_vars + artificial_lp.num_artificial_vars
        basic_vars = [i for i in range(begin, end)]
        non_basic_vars = [i for i in range(begin)]
        has_solution, solution, xb, basic_vars, non_basic_vars = _simplex(artificial_lp, basic_vars, non_basic_vars)
        
        if not has_solution:
            return has_solution, solution, xb, basic_vars, non_basic_vars        
    else:
        begin = lp.num_real_vars
        end = begin + lp.num_slack_vars
        basic_vars = [i for i in range(begin, end)]
        non_basic_vars = [i for i in range(begin)]
    
    return _simplex(lp, basic_vars, non_basic_vars)

def _simplex(lp, basic_vars, non_basic_vars):
    solved = False
    has_solution = True
    xb = None
    
    if lp.is_artificial:
        basic_vars_cp = deepcopy(basic_vars)
    
    while not solved:
        B = [[lp.a[i][j] for j in basic_vars] for i in range(len(lp.a))]
        cb = [lp.c[i] for i in basic_vars]
        xb = gauss(deepcopy(B), deepcopy(lp.b))
        lambda_t = gauss(transpose(B), deepcopy(cb))
        relative_costs = [lp.c[i] - multiply(lambda_t, column(lp.a, i)) for i in non_basic_vars]
        min_cost = argmin(relative_costs)

        if relative_costs[min_cost] >= 0.0:
            if lp.is_artificial:
                has_solution = False                        
            solved = True
        else:
            y = gauss(deepcopy(B), column(lp.a, non_basic_vars[min_cost]))            
            div = [(i, xb[i] / y[i]) for i in range(len(xb)) if y[i] > 0.0]
                
            if div:
                epsilon = min(div, key=lambda value : value[1])
                index = epsilon[0]
                basic_vars[index], non_basic_vars[min_cost] = non_basic_vars[min_cost], basic_vars[index]
                
                if lp.is_artificial:
                    a = [i for i in basic_vars if i in basic_vars_cp]
                    
                    if not a:
                        solved = True
                        non_basic_vars = [i for i in non_basic_vars if i < lp.num_real_vars + lp.num_slack_vars]
            else:
                solved = True
                has_solution = False

    if lp.is_artificial:
        return has_solution, 0.0, xb, basic_vars, non_basic_vars
    
    solution = sum([xb[i] * lp.c[basic_vars[i]] for i in range(len(basic_vars))])
    
    if lp.is_max:
        solution = -solution
    
    return has_solution, solution, xb, basic_vars, non_basic_vars