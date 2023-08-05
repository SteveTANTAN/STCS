import re
import argparse

def compute_averages_and_sums(filename):
    with open(filename, 'r') as file:
        data = file.read()
    
    attempt_pattern = "Attempt: ([\\d\\.]+)"
    attempt_matches = re.findall(attempt_pattern, data)
    total_attempts = sum(int(value) for value in attempt_matches)

    if total_attempts == 0:
        print("No attempts found, cannot compute averages")
        return
    
    fields = ["query_data_No", "time", "diameter", "total_size_of_truss",
              "total_unbalance_num", "total_percentage", "total_density"]

    sums = {}
    averages = {}
    for field in fields:
        pattern = f"{field}: ([\\d\\.]+)"
        matches = re.findall(pattern, data)
        field_sum = sum(float(value) for value in matches)
        average = field_sum / total_attempts
        sums[field] = field_sum
        averages[field] = average

    print(f"The total for Attempt is {total_attempts}")
    for field in fields:
        print(f"The total for {field} is {sums[field]}")

    print("===============================")


    print(f"The average for Attempt is 1")    
    for field in fields:
        if field == "time":
            print(f"The average for {field} is {sums[field]/100}")
        else: 
            print(f"The average for {field} is {averages[field]}")

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('filename', type=str, help='The name of the file to process')

    args = parser.parse_args()
    compute_averages_and_sums(args.filename)

if __name__ == "__main__":
    main()
