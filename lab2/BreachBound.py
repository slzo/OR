from ortools.linear_solver import pywraplp
from copy import deepcopy
import queue
from queue import PriorityQueue
import numpy as np
from fractions import Fraction
Que = queue.PriorityQueue() #优先队列
SOL = list()#最优解取值
class IP_Model(object): #求解整数规划
    def __init__(self, c, A, b, min_, max_):
        """
        c:目标函数系数
        A:系数矩阵
        b:常数项
        min_:记录x1~xn的最小值,初始时全为0
        max_:记录x1~xn的最大值,初始时全为无穷大
        """
        self._c = c
        self._A = A
        self._b = b
        self._min = min_ #下界
        self._max = max_ #上界
        self._sol = None  # 决策变量的值
        self._val = -np.inf  # 目标函数值,初始化为最小值

    def _Init(self): #初始化,将原始IP问题转化为LP问题,并加入优先队列
        v = LP_Model(self._c, self._A, self._b, self._min, self._max)
        v._solve()
        Que.put(v);

    def _solve(self):
        while (not Que.empty()): #不断选取目标函数值最小的分支问题进行处理,直到得到最优解
            lp = Que.get()
            if lp._solve()==pywraplp.Solver.INFEASIBLE : #无可行解
                continue
            #----------------有可行解-----------------
            if(lp.get_val() < self._val): #目标函数值小于下界,直接剪枝
                continue

            #判断解是否全是整数,若含非整数,继续进行分支
            sol=lp.get_sol()
            flag=0
            for i in range(len(sol)):
                if(sol[i] != int(sol[i])): #i的小数部分不为0,添加约束条件,并将2个子LP问题入优先队列
                    ch1 = deepcopy(lp)
                    ch1._max[i] = int(sol[i])
                    ch1._solve()
                    Que.put(ch1)
                    ch2 = deepcopy(lp)
                    ch2._min[i] = int(sol[i])+1
                    ch2._solve()
                    Que.put(ch2)
                    flag=1
                    break
            if(flag): #该问题进行了分支
                continue
            else: #当前问题得到了全整数解
                if(lp.get_val() == self._val):
                    SOL.append(lp.get_sol())
                    continue
                else:
                    SOL.clear();
                    self._val = lp.get_val() #更新目标函数值的下界
                    SOL.append(lp.get_sol())
    def get_sol(self): #返回可行解
        return self._sol

    def get_val(self): #返回目标函数的结果
        return self._val


class LP_Model(object): #线性规划问题求解
    def __init__(self, c, A, b, min_, max_): #各变量含义和上述IP问题相同
        self._c = c
        self._A = A
        self._b = b
        self._min = min_
        self._max = max_
        self._sol = None
        self._val = np.inf

    def _creat(self): #构造LP问题
        LP = pywraplp.Solver.CreateSolver('GLOP')
        m = len(self._b)
        n = len(self._c)
        x = [LP.NumVar(self._min[j], self._max[j], 'x[%d]' % j) for j in range(n)]
        for i in range(m):
            ct = LP.Constraint(-LP.infinity(), self._b[i])
            for j in range(n):
                ct.SetCoefficient(x[j], self._A[i][j])
        obj = LP.Objective()
        for j in range(n):
            obj.SetCoefficient(x[j], self._c[j])
        obj.SetMaximization()
        return LP, x, obj

    def _solve(self): #求解
        LP, x, obj = self._creat()
        case = LP.Solve()
        if (case == LP.OPTIMAL) | (case == LP.FEASIBLE) : #case=OPTIMAL唯一解 case=FEASIBLE无穷多解 case=INFEASIBLE无可行解
            self._sol = [x[i].solution_value() for i in range(len(x))]
            self._val = obj.Value()
        return case

    def get_sol(self): #返回可行解
        return self._sol

    def get_val(self): #返回目标函数的结果
        return self._val
    #重新定义比较函数,用于优先队列排序
    def __lt__(self, other):
        return self.get_val() > other.get_val()
    def __gt__(self, other):
        return self.get_val() < other.get_val()
    def __eq__(self, other):
        return self.get_val() == other.get_val()

def start(): #由此输入数据,并调用IP进行解决
    max_min = int(input("输入目标函数的类型,0代表max 1代表min:"))
    s = input("输入目标函数的系数(空格隔开):")
    c = np.array([Fraction(n) for n in s.split()]).astype(float)
    if(max_min):
        c*=-1
    row = int(input("输入约束条件个数:"))
    col = int(input("输入约束变量个数:"))
    A = np.ones((row,col))
    b = np.ones(row)
    print("输入约束条件矩阵(系数1~系数n,常数项,符号项(1 mean <=; 2 mean >=;0 mean = ):")
    for i in range(row):
        s = input()
        arr = [Fraction(n) for n in s.split()]
        A[i] = arr[0:col]
        opt = arr[col+1]
        b[i] = arr[col]
        if(opt == 2):
            A[i]*=-1
            b[i]*=-1
    min_ = np.zeros(col)
    max_ = np.array([np.inf for i in range(col)])
    ip = IP_Model(c,A,b,min_,max_)
    ip._Init()
    ip._solve()
    print(SOL)
    if(max_min):
        print(-1*ip.get_val())
    else:
        print(ip.get_val())



if __name__ == '__main__':
    start();
