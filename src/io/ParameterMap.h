#ifndef PARAMETERMAP_H
#define PARAMETERMAP_H

#include <cstdint>
#include <cstdio>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <type_traits>

#include "../utils/error_handling.h"

// stuff inside this namespace is only meant to be used to implement ParameterMap
namespace param_details {

/* defining a construct like this is a common workaround used to raise a compile-time error in the
 * else-branch of a constexpr-if statement. This is used to implement ``ParameterMap::try_get_``
 */
template<class> inline constexpr bool dummy_false_v_ = false;

/* Kinds of errors from converting parameters to types
 */
enum class type_err { none, generic, boolean, out_of_range };

[[noreturn]] inline void formatted_type_err_(const std::string& param, const std::string& str,
                                             const std::string& dtype, type_err type_convert_err) {
  std::string r = "";
  switch (type_convert_err) {
    case type_err::none:         r = ""; break;  // this shouldn't happen
    case type_err::generic:      r = "invalid value"; break;
    case type_err::boolean:      r = "boolean values must be \"true\" or \"false\""; break;
    case type_err::out_of_range: r = "out of range"; break;
  }
  CHOLLA_ERROR("error interpretting \"%s\", the value of the \"%s\" parameter, as a %s: %s",
               str.c_str(), param.c_str(), dtype.c_str(), r.c_str());
}

/* @{
 * helper functions that try to interpret a string as a given type.
 *
 * This returns the associated value if it has the specified type. If ``type_mismatch_is_err`` is
 * true, then the program aborts with an error if the string is the wrong type. When
 * ``type_mismatch_is_err``, this simply returns an empty result.
 */
std::optional<std::int64_t> try_int64_(const std::string& str, param_details::type_err& err);
std::optional<double> try_double_(const std::string& str, param_details::type_err& err);
std::optional<bool> try_bool_(const std::string& str, param_details::type_err& err);
std::optional<std::string> try_string_(const std::string& str, param_details::type_err& err);
/* @} */
}

/*!
 * \brief A class that provides map-like access to parameter files.
 *
 * After construction, the collection of parameters and associated values can not be mutated.
 * However, the class is not entirely immutable; internally it tracks whether parameters have been
 * accessed.
 *
 * In contrast to formats like TOML, JSON, & YAML, the parameter files don't have syntactic typing
 * (i.e. where the syntax determines formatting). In this sense, the format is more like ini files.
 * As a consequence, we internally store the parameters as strings. When we access them we need
 * to explicitly convert them to the specified type.
 */
class ParameterMap {

  struct ParamEntry { std::string param_str; bool accessed; };

private: // attributes
  std::map<std::string, ParamEntry> entries_;

public:  // interface methods

  /* Reads parameters from a parameter file and arguments.
   * 
   * \note
   * We pass in a ``std::FILE`` object rather than a filename-string because that makes testing
   * easier.
   */
  ParameterMap(std::FILE* f, int argc, char **argv);

  /* queries the number of parameters (mostly for testing purposes) */
  std::size_t size() {
    return entries_.size();
  }

  /* queries whether the parameter exists. */
  bool has_param(const std::string& param) {
    return entries_.find(param) != entries_.end();
  }

  /* queries whether the parameter exists and if it has the specified type.
   *
   * \note
   * As lThe result is always ``true``, when ``T`` is ``std::string``.
   */
  template<typename T>
  bool param_has_type(const std::string& param) {
    return try_get_<T>(param, true).has_value();
  }

  /* Retrieves the value associated with the specified parameter. If the
   * parameter does not exist or does not have the specified type, then the
   * program aborts with an error.
   *
   * \tparam The expected type of the parameter-value
   *
   * \note The name follows conventions of std::optional
   */
  template<typename T>
  T value(const std::string& param) {
    std::optional<T> result = try_get_<T>(param, false);
    if (not result.has_value()) {
      CHOLLA_ERROR("The \"%s\" parameter was not specified.", param.c_str());
    }
    return result.value();
  }

  /* @{
   * If the specified parameter exists, retrieve the associated value, otherwise return default_val. 
   * If the associated value does not have the specified type, the program aborts with an error.
   *
   * \param param The name of the parameter being queried.
   * \param default_val The value to return in case the parameter was not defined.
   *
   * \note
   * This is named after std::optional::value_or. Since the return type is commonly inferred from
   * default_val, use std::remove_cv_t to ensure nice behavior (if default_val is const).
   */
  bool value_or(const std::string& param, bool default_val) {
    return try_get_<bool>(param, false).value_or(default_val);
  }

  std::int64_t value_or(const std::string& param, int default_val) {
    return try_get_<std::int64_t>(param, false).value_or(default_val);
  }

  std::int64_t value_or(const std::string& param, long default_val) {
    return try_get_<std::int64_t>(param, false).value_or(default_val);
  }

  std::int64_t value_or(const std::string& param, long long default_val) {
    return try_get_<std::int64_t>(param, false).value_or(default_val);
  }

  double value_or(const std::string& param, double default_val) {
    return try_get_<double>(param, false).value_or(default_val);
  }

  std::string value_or(const std::string& param, const std::string& default_val) {
    return try_get_<std::string>(param, false).value_or(default_val);
  }

  std::string value_or(const std::string& param, const char* default_val) {
    return try_get_<std::string>(param, false).value_or(default_val);
  }
  /* @} */

  /* Warns about parameters that have not been accessed with the ``value`` OR ``value_or`` methods.
   *
   * \param ignore_params a set of parameter names that should never be reported as unused
   * \param abort_on_warning when true, the warning is reported as error that causes the program to
   *    abort. Default is false.
   * \param suppress_warning_msg when true, the warning isn't actually printed (this only exists for
   *    testing purposes)
   * \returns the number of unused parameters
   */
  int warn_unused_parameters(const std::set<std::string>& ignore_params,
                             bool abort_on_warning = false,
                             bool suppress_warning_msg = false) const;

private:  // private helper methods

  /* helper function template that tries to retrieve values associated with a given parameter.
   *
   * This returns the associated value if it exists and has the specified type. The returned
   * value is empty if the parameter doesn't exist. If the It can also be empty when type_abort is 
   * ``true`` and the specified type doesn't match the parameter (and is a type a parameter can
   * have).
   */
  template<typename T>
  std::optional<T> try_get_(const std::string& param, bool is_type_check) {
    auto keyvalue_pair = entries_.find(param);
    if (keyvalue_pair == entries_.end()) {
      return {};
    }

    const std::string& str = (keyvalue_pair->second).param_str;  // string associate with param

    // convert the string to the specified type and store it in out
    std::optional<T> out; // default constructed
    param_details::type_err err{}; // reports errors
    const char* dtype_name;  // for formatting errors (we use a const char* rather than a
                             // std::string so we can hold string-literals)
    if constexpr (std::is_same_v<T, bool>) {
      out = param_details::try_bool_(str, err);
      dtype_name = "bool";
    } else if constexpr (std::is_same_v<T, std::int64_t>) {
      out = param_details::try_int64_(str, err);
      dtype_name = "int64_t";
    } else if constexpr (std::is_same_v<T, double>) {
      out = param_details::try_double_(str, err);
      dtype_name = "double";
    } else if constexpr (std::is_same_v<T, std::string>) {
      out = param_details::try_string_(str, err);
      dtype_name = "string";
    } else {
      static_assert(param_details::dummy_false_v_<T>,
                    "The template type can only be bool, std::int64_t, double, or std::string.");
    }

    if (is_type_check) {
      return out; // already empty if there was a type mismatch (NEVER record parameter-access)
    } else if (err == param_details::type_err::none) {
      (keyvalue_pair->second).accessed = true; // record that we accessed the parameter
      return out;
    } else {
      param_details::formatted_type_err_(param, str, dtype_name, err);
    }
  }

};

#endif /* PARAMETERMAP_H */