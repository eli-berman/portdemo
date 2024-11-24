import java.util.Random;
import Plotter.Plotter;

public class Sorting {

	final static int SELECT_VS_QUICK_LENGTH = 12;
	final static int MERGE_VS_QUICK_LENGTH = 15;
	final static int COUNTING_VS_QUICK_LENGTH = 16;
	final static int BUBBLE_VS_MERGE_LENGTH = 12;
	final static int MERGE_VS_QUICK_SORTED_LENGTH = 11;
	final static double T = 600.0;

	/**
	 * Sorts a given array using the quick sort algorithm.
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * <p>
	 * Should run in average complexity of O(nlog(n)), and worst case complexity of O(n^2)
	 *
	 * @param arr - the array to be sorted
	 */
	public static void quickSort(double[] arr) {
		int left = 0;
		int right = arr.length - 1;
		quickSort(arr, left, right);
	}

	public static void quickSort(double[] arr, int left, int right) {

		if (left < right) {
			int partitionIndex = partition(arr, left, right);
			quickSort(arr, left, partitionIndex - 1); //before partitionIndex
			quickSort(arr, partitionIndex + 1, right); //after partitionIndex
		}
	}

	public static int partition(double[] arr, int left, int right) {
		double pivot = arr[right];
		int partitionIndex = left;

		for (int i = left; i < right; i++) {
			if (arr[i] <= pivot) {

				//swap elements of less value to the left of pivot location: arr[i] <-> arr[partitionIndex]
				double temp = arr[i];
				arr[i] = arr[partitionIndex];
				arr[partitionIndex] = temp;

				partitionIndex++;
			}
		}
		//swap pivot to final location: arr[partitionIndex] <-> pivot
		double temp = arr[right];
		arr[right] = arr[partitionIndex];
		arr[partitionIndex] = temp;

		return partitionIndex;
	}


	/**
	 * Given an array arr and an index i returns the the i'th order statstics in arr.
	 * In other words, it returns the element with rank i in the array arr.
	 * <p>
	 * At each stage the pivot is chosen to be the rightmost element of the subarray.
	 * <p>
	 * Should run in average complexity of O(n), and worst case complexity of O(n^2)
	 **/
	public static double QuickSelect(double[] arr, int i) {
		int low = 0;
		int high = arr.length - 1;
		return ithSmallest(arr, low, high, i);
	}
	public static double ithSmallest(double[] arr, int low, int high, int i){

		int partition = partition(arr, low, high);

		if (partition == i - 1) {
			return  arr[partition];
		}
		else if (partition < i - 1) {
			return ithSmallest(arr, partition + 1, high, i);
		}
		else {
			return ithSmallest(arr, low, partition - 1, i);
		}
	}


	/**
	 * Sorts a given array using the merge sort algorithm.
	 * <p>
	 * Should run in complexity O(nlog(n)) in the worst case.
	 *
	 * @param arr - the array to be sorted
	 */
	public static void mergeSort(double[] arr) {

		mergeSort(arr, 0, arr.length - 1);
	}

	public static void mergeSort(double[] arr, int left, int right) {

		if (left < right) {

			//find the middle point
			int mid = left + ((right - left) / 2);

			//sort first and second halves
			mergeSort(arr, left, mid);
			mergeSort(arr, mid + 1, right);

			//merge the sorted halves
			merge(arr, left, mid, right);
		}
	}

	public static void merge(double[] arr, int left, int mid, int right) {

		//find sizes of the two subarrays which will be merged
		int size1 = (mid - left) + 1;
		int size2 = right - mid;

		//create temporary arrays
		double[] leftArr = new double[size1];
		double[] rightArr = new double[size2];

		//copy data to temporary arrays
		for (int i = 0; i < size1; i++) {
			leftArr[i] = arr[left + i];
		}

		for (int j = 0; j < size2; j++) {
			rightArr[j] = arr[mid + 1 + j];
		}

		//merge temporary arrays
		//initial indexes of first and second subarrays
		int i = 0;
		int j = 0;

		//initial index of merged subarray array
		int k = left;
		while(i < size1 && j < size2){
			if (leftArr[i] <= rightArr[j]) {
				arr[k] = leftArr[i];
				i++;
			}
			else {
				arr[k] = rightArr[j];
				j++;
			}
			k++;
		}

		//copy remaining elements of leftArray[] if any
		while (i < size1) {
			arr[k] = leftArr[i];
			i++;
			k++;
		}

		//copy remaining elements of rightArray[] if any
		while (j < size2) {
			arr[k] = rightArr[j];
			j++;
			k++;
		}
	}



	/**
	 * Sorts a given array using bubble sort.
	 * <p>
	 * The algorithm should run in complexity O(n^2).
	 *
	 * @param arr - the array to be sorted
	 */
	public static void bubbleSort(double[] arr) {
		for (int i = 0; i < arr.length - 1; i++) {
			for (int j = 0; j < arr.length - i - 1; j++) {
				if (arr[j] > arr[j + 1]) {
					//swap next and current
					double temp = arr[j];
					arr[j] = arr[j + 1];
					arr[j + 1] = temp;

				}
			}
		}
	}

	/**
	 * Sorts a given array, using the counting sort algorithm.
	 * You may assume that all elements in the array are between 0 and k (not including k).
	 * <p>
	 * Should run in complexity O(n + k) in the worst case.
	 *
	 * @param arr - an array with positive integers
	 * @param k   - an upper bound for the values of all elements in the array.
	 */
	public static void countingSort(int[] arr, int k) {

		int n = arr.length;
		int[] output = new int[n];

		//create a count array to store number of individual digits
		//initialize array to 0
		int[] count = new int[k];
		for (int i = 0; i < k; i++) {
			count[i] = 0;
		}

		//store number count of each digit
		for (int i = 0; i < n; i++) {
			count[arr[i]]++;
		}

		//change count[i] so its position now contains true position of this digit in output array
		for (int i = 1; i <= k - 1; i++) {
			count[i] += count[i-1];
		}

		//build the output array
		//make stable by operating in reverse order
		for (int i = n - 1; i >= 0; i--) {
			output[count[arr[i]] - 1] = arr[i];
			count[arr[i]]--;
		}

		//copy output array to arr, so that arr contains sorted digits
		for (int i = 0; i < n; i++) {
			arr[i] = output[i];
		}
	}


	public static void main(String[] args) {

		countingVsQuick();
		mergeVsQuick();
		mergeVsQuickOnSortedArray();
		mergeVsBubble();
		QuickSelectVsQuickSort();

	}


	private static void countingVsQuick() {
		double[] quickTimes = new double[COUNTING_VS_QUICK_LENGTH];
		double[] countingTimes = new double[COUNTING_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < COUNTING_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumCounting = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				int[] b = new int[size];
				for (int j = 0; j < a.length; j++) {
					b[j] = r.nextInt(size);
					a[j] = b[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				countingSort(b, size);
				endTime = System.currentTimeMillis();
				sumCounting += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			countingTimes[i] = sumCounting/T;
		}
		Plotter.plot("Counting sort on arrays with elements < n", countingTimes, "Quick sort on arrays with elements < n", quickTimes);

	}

	/**
	 * Compares the merge sort algorithm against quick sort on random arrays
	 */
	public static void mergeVsQuick(){
		double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on random array", quickTimes, "merge sort on random array", mergeTimes);
	}

	/**
	 * Compares the merge sort algorithm against quick sort on pre-sorted arrays
	 */
	public static void mergeVsQuickOnSortedArray(){
		double[] quickTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		double[] mergeTimes = new double[MERGE_VS_QUICK_SORTED_LENGTH];
		long startTime, endTime;
		for (int i = 0; i < MERGE_VS_QUICK_SORTED_LENGTH; i++) {
			long sumQuick = 0;
			long sumMerge = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = j;
					b[j] = j;
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQuick += endTime - startTime;
				startTime = System.currentTimeMillis();
				mergeSort(b);
				endTime = System.currentTimeMillis();
				sumMerge  += endTime - startTime;
			}
			quickTimes[i] = sumQuick/T;
			mergeTimes[i] = sumMerge/T;
		}
		Plotter.plot("quick sort on sorted array", quickTimes, "merge sort on sorted array", mergeTimes);
	}

	/**
	 * Compares merge sort and bubble sort on random arrays
	 */
	public static void mergeVsBubble(){
		double[] mergeTimes = new double[BUBBLE_VS_MERGE_LENGTH];
		double[] bubbleTimes = new double[BUBBLE_VS_MERGE_LENGTH];
		long startTime, endTime;
		Random r = new Random();
		for (int i = 0; i < BUBBLE_VS_MERGE_LENGTH; i++) {
			long sumMerge = 0;
			long sumBubble = 0;
			for(int k = 0; k < T; k++){
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				mergeSort(a);
				endTime = System.currentTimeMillis();
				sumMerge += endTime - startTime;
				startTime = System.currentTimeMillis();
				bubbleSort(b);
				endTime = System.currentTimeMillis();
				sumBubble += endTime - startTime;
			}
			mergeTimes[i] = sumMerge/T;
			bubbleTimes[i] = sumBubble/T;
		}
		Plotter.plot("merge sort on random array", mergeTimes, "bubble sort on random array", bubbleTimes);
	}


	/**
	 * Compares the quick select algorithm with a random rank, and the quick sort algorithm.
	 */
	public static void QuickSelectVsQuickSort(){
		double[] QsortTimes = new double[SELECT_VS_QUICK_LENGTH];
		double[] QselectTimes = new double[SELECT_VS_QUICK_LENGTH];
		Random r = new Random();
		long startTime, endTime;
		for (int i = 0; i < SELECT_VS_QUICK_LENGTH; i++) {
			long sumQsort = 0;
			long sumQselect = 0;
			for (int k = 0; k < T; k++) {
				int size = (int)Math.pow(2, i);
				double[] a = new double[size];
				double[] b = new double[size];
				for (int j = 0; j < a.length; j++) {
					a[j] = r.nextGaussian() * 5000;
					b[j] = a[j];
				}
				startTime = System.currentTimeMillis();
				quickSort(a);
				endTime = System.currentTimeMillis();
				sumQsort += endTime - startTime;
				startTime = System.currentTimeMillis();
				QuickSelect(b, r.nextInt(size)+1);
				endTime = System.currentTimeMillis();
				sumQselect += endTime - startTime;
			}
			QsortTimes[i] = sumQsort/T;
			QselectTimes[i] = sumQselect/T;
		}
		Plotter.plot("quick sort with an arbitrary pivot", QsortTimes, "quick select with an arbitrary pivot, and a random rank", QselectTimes);
	}

}