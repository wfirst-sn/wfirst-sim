def read_csv(csv):
    f = open(csv)
    lines = f.read().split('\n')
    f.close()

    data = {}
    headers = lines[0].split(",")
    for item in headers:
        data[item] = []

    for line in lines[1:]:
        parsed = line.split(",")

        if len(parsed) > 2:
            assert len(parsed) == len(headers), line
            
            for i in range(len(headers)):
                item = parsed[i]
                
                try:
                    item = float(item)
                except:
                    pass

                data[headers[i]].append(item)
    return data
