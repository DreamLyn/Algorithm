package com.locaris.algorithm.tools;

public class CommonFuns {
	public static double getMax(double[] arr){
		double max = 0;
		for(int i=0;i<arr.length;i++){
			if(arr[i]>max){
				max=arr[i];
			}
		}
		return max;
	}
	public static double getMax(Double[] arr){
		double max = 0;
		for(int i=0;i<arr.length;i++){
			if(arr[i]>max){
				max=arr[i];
			}
		}
		return max;
	}
	public static double getMin(double[] arr){
		double min = arr[0];
		for(int i=1;i<arr.length;i++){
			if(arr[i]<min){
				min=arr[i];
			}
		}
		return min;
	}
	
	public static double getMin(Double[] arr){
		double min = arr[0];
		for(int i=1;i<arr.length;i++){
			if(arr[i]<min){
				min=arr[i];
			}
		}
		return min;
	}
}
