
extern_dict = {}
extern_dict["uni_adapter_count"] = 0

def worker_count_reads(seq,qual,header):
    global extern_dict
    seq = str(seq)

    universal_adapter = "AGATCGGAAGAG"
    universal_adapter = "AG"

    if universal_adapter in seq:
        extern_dict["uni_adapter_count"] = extern_dict["uni_adapter_count"] + 1


def report_counts_adapter(report,data_pack={}):
    global extern_dict

    s = "Reads with Universal Adapter" + str(extern_dict["uni_adapter_count"])
    print(s)
    report.write(s)

