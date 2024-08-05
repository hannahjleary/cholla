#include "../io/ParameterMap.h"

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <set>
#include <string>
#include <string_view>
#include <vector>

#include "../global/global.h"  // MAXLEN
#include "../io/io.h"          // chprintf
#include "../utils/error_handling.h"

[[noreturn]] void param_details::Report_TypeErr_(const std::string& param, const std::string& str,
                                                 const std::string& dtype, param_details::TypeErr type_convert_err)
{
  std::string r;
  using param_details::TypeErr;
  switch (type_convert_err) {
    case TypeErr::none:
      r = "";
      break;  // this shouldn't happen
    case TypeErr::generic:
      r = "invalid value";
      break;
    case TypeErr::boolean:
      r = R"(boolean values must be "true" or "false")";
      break;
    case TypeErr::out_of_range:
      r = "out of range";
      break;
  }
  CHOLLA_ERROR("error interpretting \"%s\", the value of the \"%s\" parameter, as a %s: %s", str.c_str(), param.c_str(),
               dtype.c_str(), r.c_str());
}

param_details::TypeErr param_details::try_bool_(const std::string& str, bool& val)
{
  if (str == "true") {
    val = true;
  } else if (str == "false") {
    val = false;
  } else {
    return param_details::TypeErr::boolean;
  }
  return param_details::TypeErr::none;
}

param_details::TypeErr param_details::try_int64_(const std::string& str, std::int64_t& val)
{
  char* ptr_end{};
  errno         = 0;  // reset errno to 0 (prior library calls could have set it to an arbitrary value)
  long long tmp = std::strtoll(str.data(), &ptr_end, 10);  // the last arg specifies base-10

  if (errno == ERANGE) {  // deal with errno first, so we don't accidentally overwrite it
    // - non-zero vals other than ERANGE are implementation-defined (plus, the info is redundant)
    return param_details::TypeErr::out_of_range;
  } else if ((str.data() + str.size()) != ptr_end) {
    // when str.data() == ptr_end, then no conversion was performed.
    // when (str.data() + str.size()) != ptr_end, str could hold a float or look like "123abc"
    return param_details::TypeErr::generic;
#if (LLONG_MIN != INT64_MIN) || (LLONG_MAX != INT64_MAX)
  } else if ((tmp < INT64_MIN) and (tmp > INT64_MAX)) {
    return param_details::TypeErr::out_of_range;
#endif
  }
  val = std::int64_t(tmp);
  return param_details::TypeErr::none;
}

param_details::TypeErr param_details::try_double_(const std::string& str, double& val)
{
  char* ptr_end{};
  errno = 0;  // reset errno to 0 (prior library calls could have set it to an arbitrary value)
  val   = std::strtod(str.data(), &ptr_end);

  if (errno == ERANGE) {  // deal with errno first, so we don't accidentally overwrite it
    // - non-zero vals other than ERANGE are implementation-defined (plus, the info is redundant)
    return param_details::TypeErr::out_of_range;
  } else if ((str.data() + str.size()) != ptr_end) {
    // when str.data() == ptr_end, then no conversion was performed.
    // when (str.data() + str.size()) != ptr_end, str could look like "123abc"
    return param_details::TypeErr::generic;
  }
  return param_details::TypeErr::none;
}

param_details::TypeErr param_details::try_string_(const std::string& str, std::string& val)
{
  // mostly just exists for consistency (every parameter can be considered a string)
  // note: we may want to consider removing surrounding quotation marks in the future
  val = str;  // we make a copy for the sake of consistency
  return param_details::TypeErr::none;
}

namespace
{  // stuff inside an anonymous namespace is local to this file

/*! Helper class that specifes the parts of a string correspond to the key and the value */
struct KeyValueViews {
  std::string_view key;
  std::string_view value;
};

/*! \brief Try to extract the parts of nul-terminated c-string that refers to a parameter name
 *  and a parameter value. If there are any issues, views will be empty optional is returned. */
KeyValueViews Try_Extract_Key_Value_View(const char* buffer)
{
  // create a view that wraps the full buffer (there aren't any allocations)
  std::string_view full_view(buffer);

  // we explicitly mimic the old behavior

  // find the position of the equal sign
  std::size_t pos = full_view.find('=');

  // handle the edge-cases (where we can't parse a key-value pair)
  if ((pos == 0) or                       // '=' sign is the first character
      ((pos + 1) == full_view.size()) or  // '=' sign is the last character
      (pos == std::string_view::npos)) {  // there is no '=' sign
    return {std::string_view(), std::string_view()};
  }
  return {full_view.substr(0, pos), full_view.substr(pos + 1)};
}

void rstrip(std::string_view& s)
{
  std::size_t cur_len = s.size();
  while ((cur_len > 0) and std::isspace(s[cur_len - 1])) {
    cur_len--;
  }
  if (cur_len < s.size()) s = s.substr(0, cur_len);
}

/*! \brief Modifies the string_view to remove trailing and leading whitespace.
 *
 *  \note
 *  Since this is a string_view, we don't actually mutate any characters
 */
void my_trim(std::string_view& s)
{
  /* Trim left side */
  std::size_t start           = 0;
  const std::size_t max_start = s.size();
  while ((start < max_start) and std::isspace(s[start])) {
    start++;
  }
  if (start > 0) s = s.substr(start);

  rstrip(s);
}

/*! This describe the current parsing status (primarily for the purpose of formatting exception messages) */
struct CurrentParseStatus {
  enum ContextKind { param_from_file, table_heading, cli_parameter};
  
  ContextKind context;
  std::string current_table_heading;
  std::string full_param_name;
  std::string_view current_param_name;

  const std::string& get_full_name() const {
    if (this->context == CurrentParseStatus::ContextKind::table_heading) return current_table_heading;
    return full_param_name;
  }

  // it's ok to be inefficient (since we are aborting)
  std::string format_err_msg(const std::string& reason) const
  {
    std::string msg = "Problem encountered while parsing the ";
    if (this->context == CurrentParseStatus::ContextKind::table_heading) {
      msg += '[';
      msg += current_table_heading;
      msg += "] parameter-table heading: "
    } else {
      msg += '"';
      msg += current_param_name;
      msg += "\" parameter ";
      if (this->context == CurrentParseStatus::ContextKind::cli_parameter) {
        msg += "from the commmand-line: "
      } else if (current_table_heading.empty()) {
        msg += "from the parameter file: "
      } else {
        msg += "under the parameter-file's [";
        msg += current_table_heading;
        msg += "] heading (aka the \"";
        msg += this->get_full_name();
        msg += "\" parameter): ";
      } 
    }
    msg += reason;
    return msg;
  }
};

/*! Helper function used to handle some parsing-related tasks that come up when considering the
 *  full name of a parameter-table or a parameter-key.
 *
 *  This does the following:
 *    1. Validates the full_name only contains allowed characters
 *    2. for a name "a.b.c.d", we step through the "a.b.c", "a.b", and "a" to
 *       (i)   confirm that no-segment is empty
 *       (ii)  ensure that the part is registered as a table
 *       (iii) ensure that the part does not collide with the name of a parameter
 *
 *  \returns An empty string if there aren't any problems. Otherwise the returned string provides an error message
 */
std::string Process_Full_Name(std::string full_name, std::set<std::string, std::less<>>& full_table_set,
                              const std::map<std::string, ParameterMap::ParamEntry>& param_entries)
{
  // first, confirm the name only holds valid characters
  std::size_t bad_value_count = 0;
  for (char ch : full_name) {
    bad_value_count +=  ((ch != '.') and (ch != '_') and (ch != '-') and not std::isalnum(ch));
  }
  if (bad_value_count > 0) {
    return "contains an unallowed character";
  }

  // now lets step through the parts of a name (delimited by the '.')
  // -> for a name "a.b.c.d", we check "a.b.c", "a.b", and "a"
  // -> specifically, we (i)   confirm that no-segment is empty
  //                     (ii)  ensure that the table is registered
  //                     (iii) ensure there aren't any collisions with parameter-names
  const std::size_t size_minus_1 = full_name.size() - 1;
  std::size_t rfind_start        = size_minus_1;
  while (true) {
    const std::size_t pos = full_name.rfind('.', rfind_start);
    if (pos == std::string_view::npos) return {};
    if (pos == size_minus_1) return "ends with a '.' character";
    if (pos == 0) return "start with a '.' character";
    if (pos == rfind_start) return "contains contiguous '.' characters";

    std::string_view table_name_prefix(full_name.data(), pos);

    // if table_name_prefix has been seen before, then we're done (its parents have been seen too)
    if (full_table_set.find(table_name_prefix) != full_table_set.end()) {
      return {};
    }

    // register table_name_prefix for the future
    std::string table_name_prefix_str(table_name_prefix);
    full_table_set.insert(table_name_prefix_str);

    if (param_entries.find(table_name_prefix_str) != param_entries.end()) {
      return "the (sub)table name collides with the existing \"" + table_name_prefix_str + "\" parameter";
    }

    rfind_start = pos - 1;
  }
}

}  // anonymous namespace

ParameterMap::ParameterMap(std::FILE* fp, int argc, char** argv)
{
  int buf;
  char *s, buff[256];

  // to provide consistent table-related behavior to TOML, we need to track the names of tables
  // (we also need to track the table-names explicitly declared in headers)
  std::set<std::string> explicit_tables;
  std::set<std::string, std::less<>> all_tables;

  CHOLLA_ASSERT(fp != nullptr, "ParameterMap was passed a nullptr rather than an actual file object");

  std::string cur_table_header{};

  /* Read next line */
  while ((s = fgets(buff, sizeof buff, fp)) != NULL) {
    /* Skip blank lines and comments */
    if (buff[0] == '\n' || buff[0] == '#' || buff[0] == ';') {
      continue;
    }

    if (buff[0] == '[') {  // here we are parsing a header like "[my_table]\n"
      std::string_view view(buff);
      rstrip(view);  // strip off trailing whitespace from the view
      if (view.back() != ']') throw std::runtime_error("problem parsing a parameter-table header");
      cur_table_header = view.substr(1, view.size() - 2);
      if (cur_table_header.size() == 0) {
        throw std::runtime_error("empty parameter-table headers (e.g. []) aren't allowed");
      }

      // confirm that we haven't seen this header before (and that there isn't a parameter with the same name)
      if (explicit_tables.find(cur_table_header) != explicit_tables.end()) {
        throw std::runtime_error("the [" + cur_table_header + "] header appears more than once");
      } else if (this->entries_.find(cur_table_header) != this->entries_.end()) {
        throw std::runtime_error("the [" + cur_table_header + "] header collides with a parameter of the same name");
      }

      std::string msg = Process_Full_Name(cur_table_header, all_tables, this->entries_);
      if (not msg.empty()) {
        throw std::runtime_error("problem encountered while parsing [" + cur_table_header + "] table header: " + msg);
      }

      // record that we've seen this headers for future checks
      explicit_tables.insert(cur_table_header);
      all_tables.insert(cur_table_header);

    } else {  // Parse name/value pair from line
      KeyValueViews kv_pair = Try_Extract_Key_Value_View(buff);
      // skip this line if there were any parsing errors (I think we probably abort with an
      // error instead, but we are currently maintaining historical behavior)
      if (kv_pair.key.empty()) continue;
      my_trim(kv_pair.value);

      if (kv_pair.key.find('.') != std::string_view::npos) {
        throw std::runtime_error(
            "the \"" + std::string(kv_pair.key) +
            "\" parameter in the contains a '.'. This isn't currently allowed in the parameter file");
      }
      std::string full_param_name = (not cur_table_header.empty()) ? (cur_table_header + '.') : std::string{};
      full_param_name += std::string(kv_pair.key);

      std::string msg = Process_Full_Name(full_param_name, all_tables, this->entries_);
      if (msg.empty()) {
        entries_[full_param_name] = {std::string(kv_pair.value), false};
      } else if (cur_table_header.empty()) {
        throw std::runtime_error("problem encountered while parsing the \"" + full_param_name + "\" parameter: " + msg);
      } else {
        throw std::runtime_error("problem encountered while parsing the \"" + std::string(kv_pair.key) +
                                 "\" parameter in the [" + cur_table_header + "] parameter-table (aka \"" +
                                 full_param_name + "\"): " + msg);
      }
    }
  }

  // Parse overriding args from command line
  for (int i = 0; i < argc; ++i) {
    // try to parse the argument
    KeyValueViews kv_pair = Try_Extract_Key_Value_View(argv[i]);
    if (kv_pair.key.empty()) continue;
    my_trim(kv_pair.value);
    std::string key_str(kv_pair.key);
    std::string msg = Process_Full_Name(key_str, all_tables, this->entries_);
    if (not msg.empty()) {
      throw std::runtime_error("problem parsing \"" + key_str + "\" parameter from the command-line: " + msg);
    }
    std::string value_str(kv_pair.value);
    chprintf("Override with %s=%s\n", key_str.c_str(), value_str.c_str());
    entries_[key_str] = {value_str, false};
  }
}

int ParameterMap::warn_unused_parameters(const std::set<std::string>& ignore_params, bool abort_on_warning,
                                         bool suppress_warning_msg) const
{
  int unused_params = 0;
  for (const auto& kv_pair : entries_) {
    const std::string& name                     = kv_pair.first;
    const ParameterMap::ParamEntry& param_entry = kv_pair.second;

    if ((not param_entry.accessed) and (ignore_params.find(name) == ignore_params.end())) {
      unused_params++;
      const std::string& value = param_entry.param_str;
      if (abort_on_warning) {
        CHOLLA_ERROR("%s/%s:  Unknown parameter/value pair!", name.c_str(), value.c_str());
      } else if (not suppress_warning_msg) {
        chprintf("WARNING: %s/%s:  Unknown parameter/value pair!\n", name.c_str(), value.c_str());
      }
    }
  }
  return unused_params;
}