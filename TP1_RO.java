
package TP;


import java.util.Scanner;
public class RO
{
	
static void GaussJordan(float A [][],int r,int s,int n,int m)
{
	float pivot =A[r][s];
    
    for (int j = 1; j <n; j++) {
    A[r][j]=A[r][j]/pivot;
    
}
    for (int i = 1; i <m; i++) {
        if(i!=s) {
        	float Ais=A[i][s];
        	for(int j=1;j<m;j++)
        	{
        		A[i][j]=A[i][j]-Ais*A[s][j];
        	}
        }
        
    }
    
}


static void getCofactor(int A[][], int temp[][], int p, int q, int n)
{
	int i = 0, j = 0;

	
	for (int r = 0; r < n; r++)
	{
		for (int s = 0; s < n; s++)
		{
			
			if (r != p && s != q)
			{
				temp[i][j++] = A[r][s];

			
				if (j == n - 1)
				{
					j = 0;
					i++;
				}
			}
		}
	}
}


static int determinant(int A[][], int n,int N)
{
	int D = 0; 

	if (n == 1)
		return A[0][0];

	int [][]temp = new int[N][N]; 

	int sign = 1; 
	for (int f = 0; f < n; f++)
	{

		getCofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1,N);


		sign = -sign;
	}

	return D;
}


static void adjoint(int A[][],int [][]adj,int N)
{
	if (N == 1)
	{
		adj[0][0] = 1;
		return;
	}

	int sign = 1;
	int [][]temp = new int[N][N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			getCofactor(A, temp, i, j, N);

			sign = ((i + j) % 2 == 0)? 1: -1;

			adj[j][i] = (sign)*(determinant(temp, N-1,N));
		}
	}
}


static boolean inverse(int A[][], float [][]inverse,int N)
{
	
	int pivot = determinant(A, N,N);
	if (pivot == 0)
	{
		System.out.print("Singular matrix, can't find its inverse");
		return false;
	}


	int [][]adj = new int[N][N];
	adjoint(A, adj,N);

	
	for (int r = 0; r < N; r++)
		for (int j = 0; j < N; j++)
			inverse[r][j] = adj[r][j]/(float)pivot;

	return true;
}


static void display(int A[][],int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			System.out.print(A[i][j]+ " ");
		System.out.println();
	}
}
static void display(float A[][],int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			System.out.printf("%.6f ",A[i][j]);
		System.out.println();
	}
}
static void equation()
{
	char []var = {'x', 'y', 'z', 'w'};
    System.out.println("Enter the number of variables in the equations: ");
    Scanner input = new Scanner(System.in);
    int n = input.nextInt();
    System.out.println("Enter the coefficients of each variable for each equations");
    System.out.println("ax + by + cz + ... = d");
    double [][]mat = new double[n][n];
    double [][]constants = new double[n][1];
   
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            mat[i][j] = input.nextDouble();
        }
        constants[i][0] = input.nextDouble();
    }
   
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            System.out.print(" "+mat[i][j]);
        }
        System.out.print("  ");
        System.out.print("  =  "+ constants[i][0]);
        System.out.println();
    }

    double inverted_mat[][] = invert(mat);
    System.out.println("The inverse is: ");
    for (int i=0; i<n; ++i)
    {
        for (int j=0; j<n; ++j)
        {
            System.out.print(inverted_mat[i][j]+"  ");
        }
        System.out.println();
    }
   
    double result[][] = new double[n][1];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            for (int k = 0; k < n; k++)
            {
                result[i][j] = result[i][j] + inverted_mat[i][k] * constants[k][j];
            }
        }
    }
    System.out.println("The product is:");
    for(int i=0; i<n; i++)
    {
        System.out.println(result[i][0] + " ");
    }
    input.close();
}
public static double[][] invert(double a[][])
{
    int n = a.length;
    double x[][] = new double[n][n];
    double b[][] = new double[n][n];
    int index[] = new int[n];
    for (int i=0; i<n; ++i)
        b[i][i] = 1;

    gaussian(a, index);

    for (int i=0; i<n-1; ++i)
        for (int j=i+1; j<n; ++j)
            for (int k=0; k<n; ++k)
                b[index[j]][k]
                        -= a[index[j]][i]*b[index[i]][k];

    for (int i=0; i<n; ++i)
    {
        x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
        for (int j=n-2; j>=0; --j)
        {
            x[j][i] = b[index[j]][i];
            for (int k=j+1; k<n; ++k)
            {
                x[j][i] -= a[index[j]][k]*x[k][i];
            }
            x[j][i] /= a[index[j]][j];
        }
    }
    return x;
}

public static void gaussian(double a[][], int index[])
{
    int n = index.length;
    double c[] = new double[n];

    for (int i=0; i<n; ++i)
        index[i] = i;

    for (int i=0; i<n; ++i)
    {
        double c1 = 0;
        for (int j=0; j<n; ++j)
        {
            double c0 = Math.abs(a[i][j]);
            if (c0 > c1) c1 = c0;
        }
        c[i] = c1;
    }

    int k = 0;
    for (int j=0; j<n-1; ++j)
    {
        double pi1 = 0;
        for (int i=j; i<n; ++i)
        {
            double pi0 = Math.abs(a[index[i]][j]);
            pi0 /= c[index[i]];
            if (pi0 > pi1)
            {
                pi1 = pi0;
                k = i;
            }
        }

        int itmp = index[j];
        index[j] = index[k];
        index[k] = itmp;
        for (int i=j+1; i<n; ++i)
        {
            double pj = a[index[i]][j]/a[index[j]][j];

            a[index[i]][j] = pj;

            for (int l=j+1; l<n; ++l)
                a[index[i]][l] -= pj*a[index[j]][l];
        }}
}

public static void main(String[] args)
{
	int A[][] = {{2,-1,0},{-1,2,-1},{0,-1,2}};


	float [][]inv = new float[3][3]; 
	System.out.print("Input matrix is :\n");
	display(A,3);


	System.out.print("\nThe Inverse is :\n");
	if (inverse(A, inv,3))
		display(inv,3);
	
	int B[][] = {{2,1},{3,1}};


	float [][]in = new float[2][2]; 
	System.out.print("Input matrix is :\n");
	display(B,2);


	System.out.print("\nThe Inverse is :\n");
	if (inverse(B, in,2))
		display(in,2);
	
	System.out.print("\n multiplication A-*A :\n");
	double c[][]=new double[3][3];  
	for(int i=0;i<3;i++){    
		for(int j=0;j<3;j++){    
		c[i][j]=0;      
		for(int k=0;k<3;k++)      
		{      
		c[i][j]+=A[i][k]*inv[k][j];      
		}  
		System.out.print(c[i][j]+" ");  
		} 
		System.out.println();   
		} 
	
	
	System.out.print("\n multiplication B-*B :\n");
	double cc[][]=new double[2][2];  
	for(int i=0;i<2;i++){    
		for(int j=0;j<2;j++){    
		cc[i][j]=0;      
		for(int k=0;k<2;k++)      
		{      
		cc[i][j]+=B[i][k]*in[k][j];      
		} 
		System.out.print(c[i][j]+" ");    
		}  
		System.out.println();    
		}   
	equation();

}
}


