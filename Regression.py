class LinearRegression:
    def __init__(self, x_count, y_count):
        self.n1 = x_count
        self.n2 = y_count
        self.num = x_count + 1
        self.X = [[1 if j == 0 else 0 for j in range(self.num)] for i in range(y_count)]
        self.Y = [[0] for _ in range(y_count)]

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
        return B

    def run(self):
        y_name = input("Input what is being estimated: ")
        print("Input parameter names:")
        x_names = [input(f"parameter {i+1}: ") for i in range(self.n1)]
        print("Input the number of observations:")
        if self.n2 <= self.n1:
            print("Number of observations should be more than number of parameters!")
            return

        for i in range(self.n2):
            for j in range(1, self.num):
                self.X[i][j] = float(input(f"{x_names[j-1]}{i+1}: "))
            self.Y[i][0] = float(input(f"{y_name}{i+1}: "))

        B = self.calculate()
        print("\nThe equation for the given data is:")
        print(f"{y_name} = ", end="")
        for i in range(len(B)):
            coef = round(B[i][0], 2)
            if i == 0:
                print(f"{coef}", end=" ")
            else:
                sign = "+" if coef >= 0 else "-"
                print(f"{sign} {abs(coef)}*{x_names[i-1]}", end=" ")
        print("\n")

        Beta = [abs(b[0]) for b in B]
        x_ordered = x_names[:]
        paired = list(zip(Beta[1:], x_ordered))
        paired.sort()
        print(f"The parameter order in deciding {y_name} is:")
        for i in range(len(paired) - 1):
            if round(paired[i][0], 2) == round(paired[i+1][0], 2):
                print(f"{paired[i][1]} = ", end="")
            else:
                print(f"{paired[i][1]} < ", end="")
        print(paired[-1][1])

        estimations = int(input("Input the number of estimations to be made: "))
        for _ in range(estimations):
            print(f"\nInput parameters to estimate {y_name}:")
            y_pred = B[0][0]
            for i in range(self.n1):
                val = float(input(f"{x_names[i]} = "))
                y_pred += B[i + 1][0] * val
            print(f"{y_name} = {round(y_pred, 2)}\n")


if __name__ == "__main__":
    print("LINEAR REGRESSION")
    n1 = int(input("Input the number of parameters: "))
    n2 = int(input("Input the number of observations: "))
    model = LinearRegression(n1, n2)
    model.run()
