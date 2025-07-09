import java.util.*;
class Linear_Regression     {
    static double[][]X,Y,transpose,product1,product2,product3,inverse,Beta;//double dimensional arrays
    static int n1,n2;
    Linear_Regression (double XY[][],int x,int y)    {//parameterized constructor
        transpose = new double[XY[0].length][XY.length];
        product1= new double[XY[0].length][XY[0].length];
        product2=new double[XY[0].length][1];
        product3= new double[XY[0].length][1];
        Beta= new double[XY.length][1];
        inverse= new double[XY[0].length][XY[0].length];
        X = new double[XY.length][XY[0].length];
        Y = new double[XY.length][1];
        n1=x;
        n2=y;
    }

    static double[][] transpose(double[][] matrix)  {//finding the transpose of the matrix
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                transpose[j][i] = matrix[i][j];//storing the transpose in array transpose
        return transpose;//returning the transpose
    }

    static double[][] multiply_matrices(double[][] X,double[][]Y)   {//finding the product of the two matrices
        double P[][]= new double[X.length][Y[0].length];
        for(int i=0;i<X.length;i++)
            for(int j=0;j<Y[0].length;j++)
                for (int c=0;c<Y.length;c++)
                    P[i][j]+=X[i][c]*Y[c][j];
        return P;//returning the product
    }

    static void getCofactor(double A[][],double temp[][],int p,int q,int n) {//finding the cofactor matrix of A
        int i = 0, j = 0;
        for (int row = 0; row < n; row++)   {
            for (int col = 0; col < n; col++)   {
                if (row != p && col != q)   {
                    temp[i][j++] = A[row][col];
                    if (j == n - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    static int determinant(double A[][], int n) {//finding the determinant 
        int D = 0;//intialising D
        if (n == 1)//base case
            return (int)A[0][0];//determinant if matrix is 1X1
        double [][]temp = new double[X[0].length][X[0].length]; 
        int sign = 1;
        for (int f = 0; f < n; f++) {
            getCofactor(A, temp, 0, f, n);//updating temp as the cofactor of A by invoking getCofactor function
            int B=(int)A[0][f];
            int AB=determinant(temp, n - 1);//recursive case
            D += sign * B *AB;//updating the determinant
            sign = -sign;
        }
        return D;//returning the determinant
    }

    static void adjoint(double A[][],double [][]adj)    {//finding the adjoint of the matrix
        if (X.length == 1)  {
            adj[0][0] = 1;
        }
        int sign = 1;
        double [][]temp = new double[X[0].length][X[0].length];
        for (int i = 0; i < X[0].length; i++)   {
            for (int j = 0; j < X[0].length; j++)   {
                getCofactor(A, temp, i, j, X[0].length);
                sign = ((i + j) % 2 == 0)? 1: -1;//alternating sign between 1 and -1
                adj[j][i] = (sign)*(determinant(temp, X[0].length-1));
            }
        }
    }

    static double[][] inverse(double A[][]) {//finding inverse of matrix
        int det = determinant(A, X[0].length);
        double [][]adj = new double[X[0].length][X[0].length];
        adjoint(A, adj);
        for (int i = 0; i < X[0].length; i++)
            for (int j = 0; j < X[0].length; j++)
                inverse[i][j] = adj[i][j]/(float)det;// inverse = adjoint matrix/determinant
        return inverse;
    }

    static double[][] calculation() {//calculating the linear regression using the formula
        //forumla : beta(coordinates of equation) = product(inverse( product( transpose(X) and X)) and product(transpose(X) and Y))
        transpose=transpose(X);
        product1=multiply_matrices(transpose,X);
        product2=multiply_matrices(transpose,Y);
        inverse=inverse(product1);
        product3=multiply_matrices(inverse,product2);
        return product3;
    }

    public static void main(String args[])  {
        System.out.println("LINEAR REGRESSION ");
        Scanner sc= new Scanner(System.in);
        System.out.print("Input what is being estimated :");//Taking inputs
        String y= sc.nextLine();
        System.out.print("Input the number of parameters :");
        int n1= sc. nextInt();
        int num=n1+1;
        String x2[]=new String[n1];
        String x[]=new String[n1];
        double x1[]=new double[num];
        double y1=0.0;
        System.out.println("Input parameter names :");
        x[0]=sc.nextLine();
        for(int i=0;i<n1;i=i+1) {
            System.out.print("parameter"+(i+1)+" : ");
            x[i]=sc.nextLine();
            x2[i]=x[i];
        }
        System.out.print("Input the number of observations :");
        int n2= sc. nextInt();
        if(n2<=n1)
            System.out.println("Number of observations should be more than number of parameters!");
        else    {
            X = new double[n2][num];
            Y = new double[n2][1];
            Linear_Regression  obj= new Linear_Regression (X,n1,n2);
            for(int i=0;i<n2;i++)
                X[i][0]=1;
            for(int i=0;i<n2;i++) 
                for(int j=1;j<=num;j++)  {
                    if(j==(num)) {
                        System.out.print(y+(i+1)+" : ");
                        Y[i][0]=sc.nextDouble();
                        System.out.println();
                    }
                    else    {
                        System.out.print(x[j-1]+(i+1)+" : ");
                        X[i][j]=sc.nextDouble();
                    }
                }
            product3=obj.calculation();//finding the coefficients of the equation by invoking calculation method
            for(int i=0;i<n1;i++)
                System.out.print(x[i]+"\t ");
            //printing the given values
            System.out.print(y);
            System.out.println();
            for(int i=0;i<n2;i++)   {
                for(int j=1;j<=num;j++)  {
                    if(j==(num))
                        System.out.print(Y[i][0]+"\t");
                    else
                        System.out.print(X[i][j]+"\t");
                }
                System.out.println("");
            }
            //printing the equation
            System.out.print("The equation for the given data is :");
            System.out.print(y+"=");
            Beta[0][0]=product3[0][0];
            if(product3[1][0]>=0)
                System.out.print(Math.round(Beta[0][0] * 100.0) / 100.0+" + ");
            else
                System.out.print(Math.round(Beta[0][0] * 100.0) / 100.0+" - ");
            for(int i=1;i<num-1;i++) {
                if(i!=0)
                    Beta[i][0]=Math.abs(product3[i][0]);
                if(product3[i+1][0]>=0)
                    System.out.print(Math.round(Beta[i][0] * 100.0) / 100.0+" * "+x[i-1]+" + ");
                else
                    System.out.print(Math.round(Beta[i][0] * 100.0) / 100.0+" * "+x[i-1]+" - ");
            }
            Beta[num-1][0]=Math.abs(product3[num-1][0]);
            System.out.print(Math.round(Beta[n1][0] * 100.0) / 100.0+" * "+x[n1-1]+" ");
            System.out.println("\n");
            //printing the order in which the given parameters affct the final value
            if(num>1)    {
                System.out.print("The parameter order in deciding "+y+" is :");
                for (int i = 1; i < num-1; i++)  {
                    int min_idx = i;
                    for (int j = i+1; j < num; j++)
                        if (Beta[j][0] < Beta[min_idx][0])
                            min_idx = j;
                    double temp = Beta[min_idx][0];
                    Beta[min_idx][0] = Beta[i][0];
                    Beta[i][0] = temp;
                    String Temp = x[min_idx-1];
                    x[min_idx-1] = x[i-1];
                    x[i-1]= Temp;
                }
                for(int i=0;i<n1-1;i++)
                    if (Math.round(Beta[i+1][0] * 100.0) / 100.0 == Math.round(Beta[i+2][0] * 100.0) / 100.0)
                        System.out.print(x[i]+" = ");
                    else
                        System.out.print(x[i]+" < ");
                System.out.println(x[n1-1]+"\n");
            }
            //taking values for new estimations
            System.out.print("Input the number of estimations to be made :");
            int ab=sc.nextInt();
            for(int j=0;j<ab;j++)   {
                System.out.println("Input parameters to estimate "+y);
                y1=product3[0][0];
                for(int i=0;i<n1;i++)   {

                    System.out.print(x2[i]+" = ");
                    x1[i]=sc.nextDouble();
                    y1+=product3[i+1][0]*x1[i];
                }
                //printing the values for given data
                System.out.println(y+" = "+Math.round(y1* 100.0) / 100.0+"\n");
            }
        }
    }
}//end of class