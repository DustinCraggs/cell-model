#include <iostream>

__global__
void test(int n, int *a, int *b, int *c) {
	int index = threadIdx.x;
	int stride = blockDim.x;
	for (int i = index; i < n; i += stride) {
		c[i] = a[i] + b[i];
	}
}

int main() {
	int n = 10000;
	int *a, *b, *c;
	cudaMallocManaged(&a, n * sizeof(int));
	cudaMallocManaged(&b, n * sizeof(int));
	cudaMallocManaged(&c, n * sizeof(int));

	for (int i = 0; i < n; i++) {
		std::cout << i << std::endl;
		a[i] = 1;
		b[i] = 2;
		c[i] = 0;
	}

	int blockSize = 256;
	int numBlocks = (n + blockSize - 1) / blockSize;
	test<<<numBlocks, blockSize>>>(n, a, b, c);

	std::cout << c[100] << std::endl;

}