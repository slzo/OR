#include <iostream>
#include <string>
using namespace std;


int op; // op = 1 = max; op = 0 = min
int row, col; //约束条件系数矩阵维数
int tet; // tet = row+col ,若标准化之后tet > col 需要2阶段
double **A, *B, *O; //约束条件系数矩阵，常量矩阵，符号矩阵
double	*C; //目标函数系数矩阵
double *Segama; //检验值
int *base; //基向量index
int stage = 1;

void Cal_Segama() { //计算检验值
	for(int i = 0; i < col; i++) {
			double tmp = C[i];
			for(int j = 0; j < row; j++) 
					tmp -= A[j][i]*C[base[j]];
			Segama[i] = tmp;
	}
	return ;
}

void Deal() { //转化为标准形式
	for(int i = 0; i< row; i++) { //常数项转正
		if(B[i] < 0) {
			for(int j = 0; j < col; j++)
				A[i][j] *= -1;
			O[i] *= -1;
		}
	}
	for(int i = 0; i < row; i++) { //添加松弛变量 全部转化为 =
			if(O[i] == 0) { // =
				continue;
			} 
			else if(O[i] == 1) { // <=
				C[col] = 0;
				A[i][col++] = 1;
				O[i] = 0;
			}
			else { // >=
				C[col] = 0;
				A[i][col++] = -1;
				O[i] = 0;
			}
	}
	if(!op) { //min -> max
		for(int i = 0; i < col; i++)
			   C[i] *= -1;
	}	
	return ;
}
	
void Getdata() {
	cin >> op; // op = 1 = max ; op = 0 = min
	cin >> row >> col;
	tet = row+col+row;

	C = new double[tet]; //目标函数系数向量
	for(int i = 0; i < col; i++)
		cin >> C[i];

	A = new double*[row];
	for(int i = 0; i < row; i++)
			A[i] = new double[tet];
	B = new double[row];
	O = new double[row];
	for(int i = 0; i < row; i++) { //输入约束条件矩阵
		for(int j = 0; j < col; j++)
			cin >> A[i][j];
		for(int j = col; j < tet; j++)
			A[i][j] = 0;
		cin >> B[i];
		cin >> O[i];
	}
	Deal();
	return ;
}

void Print() {
    cout << "------------------------------------------------------------" << endl;
	cout << "          C           |   ";
	for(int i = 0; i < col; i++)
			cout << C[i] << "  ";
    cout << endl << "------------------------------------------------------------" << endl;
    cout << "  Cb  |  Xb   |   B   |   ";
	for(int i = 0; i < col; i++)
		cout << "x" << i+1 << "  ";
    cout << endl << "------------------------------------------------------------" << endl;
	for(int i = 0; i < row; i++) {
			cout << "   " << C[base[i]] << "  |  ";
			cout << "x" << base[i]+1 << "  |  ";
			cout << B[i] << "  |  ";
			for(int j = 0; j < col; j++)
					cout << A[i][j] << "  " ;
			cout << endl;
	}
    cout << "------------------------------------------------------------" << endl;
	cout << "       Sigama       |";
	cout.width(6);
	for(int i = 0; i < col; i++)
			cout << Segama[i] << "  ";
	cout << endl;
	return ;
}

int Simplex() {
		Segama = new double[col];
		int step = 1;
		while(1) {
			Cal_Segama();
			cout << endl << "Stage" << stage << ".  Step" << step++ <<": " << endl;
			Print();
			double maxs = 0;
		   	int index = -1;
			for(int i = 0; i < col; i++)
					if(Segama[i] > maxs) {
							index = i; //换入元素的下标
							maxs = Segama[i];
					}
			if(maxs == 0) {//所有的检验值均小于0
				int counter = 0;
				for(int i = 0; i < col; i++)
					if(Segama[i] == 0)
						counter++;
				if(counter > row)
					return 3;
				return 1;
			}
			int outdex = -1;
	   		double maxb = 1000;
			for(int i = 0; i < row; i++) { //遍历寻找换出元素
				if(A[i][index] > 0) {
						double tmp = B[i]/A[i][index];
						if(tmp < maxb) {
								maxb = tmp;
								outdex = i;
						}
				}
			}
			if(outdex == -1) //有Segama>0,Pi<0 有无界解
					return 2;
			//进行基向量变换
			base[outdex] = index; //基向量换入

			//更新单纯性表
			B[outdex] /= A[outdex][index];
			for(int i = 0; i < col; i++) {
				if( i == index)
					continue;
				A[outdex][i] /= A[outdex][index];
			}
			A[outdex][index] = 1;
			for(int i = 0; i < row; i++) {
					if(i == outdex) 
							continue;
					for(int j = 0; j < col; j++) {
						   	if(j == index)
								continue;	
							A[i][j] -= A[i][index]*A[outdex][j];
					}
					B[i] -= A[i][index]*B[outdex];
					A[i][index] = 0;
			}
		}
		return 1;
}

int Base() { //获取基向量并判断是否存在单位矩阵 -> 是否需要二阶段法
	int flag = 1;
	base = new int[row];
	for(int i = 0; i < row; i++) {
		int flag1 = 0;
		for(int j = 0; j < col; j++) {
			if(A[i][j] == 1) {
				flag1 = 1;
				for(int k = 0; k < row; k++) {
					if(k == i) 
						continue;
					if(A[k][j]) {
						flag1 = 0;
						break;
					}
				}
				if(flag1) { //选为基向量
					base[i] = j;
					break;
				}
			}
		}
		if(!flag1) { //添加人工变量
				flag = 0; //需要2阶段
				A[i][col] = 1;
				base[i] = col++;
		}
	}
	return flag;
}
void ans(int re){
		double z = 0;
		if(re == 1) { //case1 有唯一最优解
			cout << "The optimal solution of the problrm is" << endl;
			cout << "x* = [";
			for(int i = 0; i < col; i++){
				int o = 0;
				for(int j = 0; j < row; j++){
					if(i == base[j]) {
						cout << B[j] << " ";
						z += B[j]*C[i];
						o = 1;
						break;
					}
				}
				if(o) continue;
				else cout << 0 << " ";
			}
			cout << "]" << endl;
			cout << "z* = ";
		    if(!op)
				cout << z*-1 << endl;
			else cout << z << endl;
		}
		else if(re == 2){ //case2：有无界解
			cout << "The optimal solution of the problem is unbounded" << endl;
			cout << "x* is unbounded" << endl;
			cout << "z* is unbounded" << endl;
		}
		if(re == 3) { //case3 ：有无限多最优解
				for(int i = 0; i < col; i++){
					for(int j = 0; j < row; j++){
						if(i == base[j]) {
							z += B[j]*C[i];
							break;
						}
				}
			}
			cout << "The number of optimal solution is unlimited" << endl;
			cout << "z* is  " ;
			if(!op)
				cout << z*-1 <<endl;
			else cout << z << endl;
		}

	return;
}
void Solve(){
	Getdata();
	int tmp = col;
	if(Base())  //直接单纯性法
		ans( Simplex() );

	else { //2阶段单纯性法
		double *C1 = new double[col]; //记录初始的参数
		for(int i = 0; i < tmp; i++) {
				C1[i] = C[i];
				C[i] = 0;
		}
		for(int i = tmp; i < col; i++) {
				C1[i] = C[i];
				C[i] = -1;
		}
		int re = Simplex();
		for(int i = tmp; i < col; i++) //检查人工变量是否全为0 不全0则无解
			for(int j = 0; j < row; j++)
				if(base[j] == i && B[j]){
					cout << "The problem doesn't has a feasible solution." << endl;
					return;
				}
		ans(re);
		col = tmp;
		for(int i = 0; i < col; i++)
				C[i] = C1[i];
		stage = 2;
		re = Simplex();
		ans(re);
	}
	return ;
}

int main() {
	Solve();
	return 0;
}
