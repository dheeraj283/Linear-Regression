import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # for 3D plots
import numpy as np


class LinearRegression:
    def __init__(self, x_count, y_count):
        self.n1 = x_count
        self.n2 = y_count
        self.num = x_count + 1
        self.X = [[1 if j == 0 else 0 for j in range(self.num)] for i in range(y_count)]
        self.Y = [[0] for _ in range(y_count)]
        self.x_names = []
        self.y_name = ""
        self.B = []

    def transpose(self, matrix):
        return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

    def multiply_matrices(self, A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

    def get_cofactor(self, A, p, q, n):
        temp = []
        for i in range(n):
            if i == p:
                continue
            row = []
            for j in range(n):
                if j == q:
                    continue
                row.append(A[i][j])
            temp.append(row)
        return temp

    def determinant(self, A):
        n = len(A)
        if n == 1:
            return A[0][0]
        D = 0
        sign = 1
        for f in range(n):
            temp = self.get_cofactor(A, 0, f, n)
            D += sign * A[0][f] * self.determinant(temp)
            sign = -sign
        return D

    def adjoint(self, A):
        n = len(A)
        if n == 1:
            return [[1]]
        adj = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                temp = self.get_cofactor(A, i, j, n)
                sign = 1 if (i + j) % 2 == 0 else -1
                adj[j][i] = sign * self.determinant(temp)
        return adj

    def inverse(self, A):
        det = self.determinant(A)
        if det == 0:
            raise ValueError("Singular matrix, can't find inverse")
        adj = self.adjoint(A)
        return [[adj[i][j] / det for j in range(len(adj))] for i in range(len(adj))]

    def calculate(self):
        XT = self.transpose(self.X)
        XTX = self.multiply_matrices(XT, self.X)
        XTY = self.multiply_matrices(XT, self.Y)
        XTX_inv = self.inverse(XTX)
        B = self.multiply_matrices(XTX_inv, XTY)
        self.B = B
        return B

    def plot(self):
        if self.n1 == 1:
            # 2D Line plot
            x_vals = [self.X[i][1] for i in range(self.n2)]
            y_vals = [self.Y[i][0] for i in range(self.n2)]
            x_range = np.linspace(min(x_vals), max(x_vals), 100)
            y_pred = self.B[0][0] + self.B[1][0] * x_range

            plt.scatter(x_vals, y_vals, color="blue", label="Original Data")
            plt.plot(x_range, y_pred, color="red", label="Regression Line")
            plt.xlabel(self.x_names[0])
            plt.ylabel(self.y_name)
            plt.legend()
            plt.title("Linear Regression Plot")
            plt.grid(True)
            plt.show()

        elif self.n1 == 2:
            # 3D Plane plot
            x_vals = [self.X[i][1] for i in range(self.n2)]
            y_vals = [self.X[i][2] for i in range(self.n2)]
            z_vals = [self.Y[i][0] for i in range(self.n2)]

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(x_vals, y_vals, z_vals, color='blue', label='Original Data')

            x_range = np.linspace(min(x_vals), max(x_vals), 20)
            y_range = np.linspace(min(y_vals), max(y_vals), 20)
            x_mesh, y_mesh = np.meshgrid(x_range, y_range)
            z_mesh = (self.B[0][0] +
                      self.B[1][0] * x_mesh +
                      self.B[2][0] * y_mesh)

            ax.plot_surface(x_mesh, y_mesh, z_mesh, alpha=0.5, color='red')
            ax.set_xlabel(self.x_names[0])
            ax.set_ylabel(self.x_names[1])
            ax.set_zlabel(self.y_name)
            plt.title("Linear Regression Surface")
            plt.legend()
            plt.show()

        else:
            print("Plotting is only supported for 1 or 2 parameters.")

    def run(self):
        self.y_name = input("Input what is being estimated: ")
        print("Input parameter names:")
        self.x_names = [input(f"parameter {i+1}: ") for i in range(self.n1)]

        if self.n2 <= self.n1:
            print("Number of observations should be more than number of parameters!")
            return

        for i in range(self.n2):
            for j in range(1, self.num):
                self.X[i][j] = float(input(f"{self.x_names[j-1]}{i+1}: "))
            self.Y[i][0] = float(input(f"{self.y_name}{i+1}: "))

        B = self.calculate()

        # Print equation
        print("\nThe equation for the given data is:")
        print(f"{self.y_name} = ", end="")
        for i in range(len(B)):
            coef = round(B[i][0], 2)
            if i == 0:
                print(f"{coef}", end=" ")
            else:
                sign = "+" if coef >= 0 else "-"
                print(f"{sign} {abs(coef)}*{self.x_names[i-1]}", end=" ")
        print("\n")

        # Plotting
        self.plot()

        # Estimation
        estimations = int(input("Input the number of estimations to be made: "))
        for _ in range(estimations):
            print(f"\nInput parameters to estimate {self.y_name}:")
            y_pred = B[0][0]
            for i in range(self.n1):
                val = float(input(f"{self.x_names[i]} = "))
                y_pred += B[i + 1][0] * val
            print(f"{self.y_name} = {round(y_pred, 2)}\n")


if __name__ == "__main__":
    print("LINEAR REGRESSION")
    n1 = int(input("Input the number of parameters (1 or 2 recommended for plotting): "))
    n2 = int(input("Input the number of observations: "))
    model = LinearRegression(n1, n2)
    model.run()
